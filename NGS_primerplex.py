#!/usr/bin/python3
# This script constructs primers for multiplex NGS panels

import argparse,re,os,sys,primer3,math,logging,pysam,xlrd
from multiprocessing.pool import ThreadPool
from Bio import SeqIO,Seq,pairwise2
import subprocess as sp
from operator import itemgetter
import xlsxwriter as xls
import networkx as nx
import networkx.algorithms.clique as clique
from itertools import islice
import myvariant
from collections import Counter
import random

global thisDir
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

# Section of functions
def readInputFile(regionsFile):
    allRegions={}
    regionsNames={}
    regionNameToChrom={}
    regionsCoords={}
    regionNameToMultiplex={}
    regionNameToPrimerType={}
    uniquePointRegions=0
    totalPointInputRegions=0
    for string in regionsFile:
        if ('Chrom' in string or
            'Start\tEnd' in string or
            string=='' or string=='\n'):
            continue            
        cols=string.replace('\n','').split('\t')
        chrom=cols[0].replace('chr','')
        if chrom=='X':
            chrom='23'
        elif chrom=='Y':
            chrom='24'
        try:
            regStart=int(cols[1])
            regEnd=int(cols[2])
        except ValueError:
            print('ERROR: Incorrect format of input file!')
            print('It should have the following format:')
            print('Chromosome{Tab}Start_Position{Tab}End_Position{Tab}Amplicon_Name{Tab}\n'
                  'Desired_Multiplex_Numbers(optional){Tab}Type_Of_Primers(only left/only right/both)(optional){Tab}'
                  'Use_Whole_Region(optional)')
            print('But your file have the following format:')
            print(string)
            exit(1)
        # Check that region is written in the right way
        if regStart>regEnd:
            temp=regStart
            regStart=regEnd
            regEnd=temp
        regionName=cols[3]
        # We split any region onto list of poses
        for i in range(regStart,regEnd+1):
            totalPointInputRegions+=1
            try:
                if int(chrom) not in regionsCoords.keys() or i not in regionsCoords[int(chrom)]:
                    if regionName not in regionsNames.keys():
                        regionsNames[regionName]=[regionName+'_1']
                    else:
                        prevNum=regionsNames[regionName][-1].split('_')[-1]
                        regionsNames[regionName].append(regionName+'_'+str(int(prevNum)+1))
                    curRegionName=regionsNames[regionName][-1]
                    regionNameToChrom[curRegionName]=chrom
                    if len(cols)>4 and cols[4]!='':
                        regionNameToMultiplex[curRegionName]=cols[4].split(',')
                    if len(cols)>5 and cols[5]!='':
                        if cols[5] not in ['L','R','B']:
                            print('ERROR! Unknown type of primers is necessary to be designed for the following line of input file:')
                            print(string)
                            print('This value can be only "L" (only left primer), "R" (only right primer), "B" (both primers) or nothing (both primers)')
                            exit(1)
                        regionNameToPrimerType[curRegionName]=cols[5]
                    if len(cols)>6 and cols[6]!='':
                        if cols[6] not in ['W']:
                            print('ERROR! Unknown value for "using whole region" for the following line of input file:')
                            print(string)
                            print('This value can be only "W" (do not split this region onto several amplicons), or nothing (region can be splitted onto several amplicons)')
                            exit(1)
                        if cols[6]=='W':
                            endShift=regEnd-regStart
                        else:
                            endShift=0
                    else:
                        endShift=0
                    if chrom not in allRegions.keys():
                        allRegions[chrom]={curRegionName:[chrom,i,i+endShift,curRegionName]}
                        regionsCoords[int(chrom)]=[i]
                        uniquePointRegions+=1
                    else:
                        # Addtionally check that all input regions are unique
                        if [chrom,i,i+endShift,curRegionName] not in allRegions[chrom].values():
                            allRegions[chrom][curRegionName]=[chrom,i,i+endShift,curRegionName]
                            regionsCoords[int(chrom)].append(i)
                            uniquePointRegions+=1
                    if endShift>0:
                        break
            except ValueError:
                print('ERROR: Incorrect format of input file!')
                print('It should have the following format:')
                print('Chromosome{Tab}Start_Position{Tab}End_Position{Tab}Amplicon_Name{Tab}\n'
                      'Desired_Multiplex_Numbers(optional){Tab}Type_Of_Primers(only left/only right/both)(optional){Tab}'
                      'Use_Whole_Region(optional)')
                print('But your file have the following format:')
                print(string)
                exit(1)
    # Sort dict with regions coords
    for chrom in regionsCoords.keys():
        regionsCoords[chrom]=sorted(regionsCoords[chrom])
    print(' # Total number of input point regions:',totalPointInputRegions)
    print(' # Number of unique input point regions:',uniquePointRegions)
    logger.info(' # Total number of input point regions: '+str(totalPointInputRegions))
    logger.info(' # Number of unique input point regions: '+str(uniquePointRegions))
    return(allRegions,regionsNames,regionsCoords,regionNameToChrom,regionNameToMultiplex,regionNameToPrimerType)

def readPrimersFile(primersFile):
    wb=xlrd.open_workbook(primersFile)
    ws=wb.sheet_by_index(0)
    internalPrimers=[]
    amplifiedRegions={}
    for i in range(ws.nrows):
        row=ws.row_values(i)
        if row[0]=='': break
        if i==0:
            # Check that input file has necessary format
            if row[:16]!='#	Left_Primer_Seq	Right_Primer_Seq	Amplicon_Name	Chrom	Amplicon_Start	Amplicon_End	Amplicon_Length	Amplified_Block_Start	Amplified_Block_End	Left_Primer_Tm	Right_Primer_Tm	Left_Primer_Length	Right_Primer_Length	Left_GC	Right_GC'.split('\t'):
                print('ERROR! Input file with primers has incorrect format. You can use only files with primers that has format of NGS-primerplex')
                logger.error('Input file with primers has incorrect format. You can use only files with primers that has format of NGS-primerplex')
                exit(1)
            continue
        if row[4] not in amplifiedRegions.keys():
            amplifiedRegions[row[4]]=set(range(int(row[8]),int(row[9])+1))
        else:
            amplifiedRegions[row[4]].update(set(range(int(row[8]),int(row[9])+1)))
        internalPrimers.append(row[1:16])
    return(internalPrimers,amplifiedRegions)

def createPrimer3_parameters(pointRegions,args,species='human',
                             designedInternalPrimers=None,
                             regionNameToPrimerType=None,
                             regionNameNeedToBeWholeLen=None,
                             amplToChrom={},
                             amplToStartCoord={}):
    templatePrimerTags={'PRIMER_TASK':'generic',
                        'PRIMER_PICK_LEFT_PRIMER':1,
                        'PRIMER_PICK_RIGHT_PRIMER':1,
                        'PRIMER_MAX_NS_ACCEPTED':0,
                        'PRIMER_SALT_CORRECTIONS':1,
                        'PRIMER_TM_FORMULA':1,
                        'PRIMER_DNA_CONC':args.primerConc, # primer, in nM
                        'PRIMER_SALT_DIVALENT':args.dvConc, # Mg2+, in mM
                        'PRIMER_SALT_MONOVALENT':args.mvConc, # K+ or NH4+, in mM
                        'PRIMER_DNTP_CONC':args.dntpConc, # sum concentration of each dNTP, in mM
                        'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT':1,
                        'PRIMER_WT_SELF_END':0.5,
                        'PRIMER_PAIR_WT_DIFF_TM':0.5,
                        'PRIMER_WT_GC_PERCENT_GT':0.5,
                        'PRIMER_WT_GC_PERCENT_LT':0.5,
                        'PRIMER_PAIR_WT_PRODUCT_SIZE_GT':2,
                        'PRIMER_PAIR_WT_PRODUCT_SIZE_LT':2,
                        'PRIMER_EXPLAIN_FLAG':1}
    primer3Params={}
    primerTags=templatePrimerTags.copy()
    # Writing primer design parameters that user defined
    primerTags['PRIMER_MIN_SIZE']=args.minPrimerLen
    primerTags['PRIMER_MAX_SIZE']=args.maxPrimerLen
    primerTags['PRIMER_OPT_SIZE']=args.optPrimerLen
    primerTags['PRIMER_MIN_TM']=args.minPrimerMelt
    primerTags['PRIMER_MAX_TM']=args.maxPrimerMelt
    primerTags['PRIMER_OPT_TM']=args.optPrimerMelt
    primerTags['PRIMER_MIN_GC']=args.minPrimerGC
    primerTags['PRIMER_MAX_GC']=args.maxPrimerGC
    primerTags['PRIMER_OPT_GC_PERCENT']=args.optPrimerGC
    primerTags['PRIMER_MAX_END_GC']=args.maxPrimerEndGC
    primerTags['PRIMER_MIN_END_GC']=args.minPrimerEndGC
    primerTags['PRIMER_MAX_POLY_X']=args.maxPrimerPolyN
    primerTags['PRIMER_MAX_SELF_END']=args.maxPrimerComplEndTh
    primerTags['PRIMER_MAX_SELF_ANY_TH']=args.maxPrimerComplAnyTh
    primerTags['PRIMER_MAX_HAIRPIN_TH']=args.maxPrimerHairpinTh
    primerTags['PRIMER_NUM_RETURN']=args.primernum1
##    primerTags['PRIMER_PAIR_MAX_COMPL_ANY_TH']=args.maxPrimerComplAnyTh ## Use of this value lead to mistake in primer3-py module. It begins to incorrectly calculate Tm
##    primerTags['PRIMER_PAIR_MAX_COMPL_END_TH']=args.maxPrimerComplEndTh ## Use of this value lead to mistake in primer3-py module. It begins to incorrectly calculate Tm
    # If we are going to design external primers for already created internal primers
    if designedInternalPrimers:
        for internalPrimers in designedInternalPrimers.values():
            leftPrimer,rightPrimer,amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,leftPrimerTm,rightPrimerTm,leftPrimerLen,rightPrimerLen,leftGC,rightGC=internalPrimers
            curRegionName=amplName[:amplName.rfind('_')]
            if regionNameToPrimerType is None:
                primerTags['PRIMER_PICK_LEFT_PRIMER']=1
                primerTags['PRIMER_PICK_RIGHT_PRIMER']=1
            elif curRegionName in regionNameToPrimerType.keys():
                if regionNameToPrimerType[curRegionName]=='L':
                    primerTags['PRIMER_PICK_LEFT_PRIMER']=1
                    primerTags['PRIMER_PICK_RIGHT_PRIMER']=0
                elif regionNameToPrimerType[curRegionName]=='R':
                    primerTags['PRIMER_PICK_LEFT_PRIMER']=0
                    primerTags['PRIMER_PICK_RIGHT_PRIMER']=1
                elif regionNameToPrimerType[curRegionName]=='B':
                    primerTags['PRIMER_PICK_LEFT_PRIMER']=1
                    primerTags['PRIMER_PICK_RIGHT_PRIMER']=1
                else:
                    print('ERROR! Unknown type of primers is necessary to be designed for the following region:')
                    print(amplName)
                    print('This value can be only "L" (only left primer), "R" (only right primer), "B" (both primers) or nothing (both primers)')
                    exit(1)
            amplToChrom[amplName]=chrom
            primer3Params[amplName]=[]
            primerTags['PRIMER_PRODUCT_SIZE_RANGE']=[[amplLen+2*args.minPrimerShift,args.maxExtAmplLen]]
            primerTags['PRIMER_PRODUCT_OPT_SIZE']=args.optExtAmplLen
            # We can shift our external amplicon on the (maximal length of extAmpl) - (amplLen)-(minPrimerShift)
            chrTargetSeqStart=amplBlockEnd-(args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen)
            # In COSMIC database chromosome X and Y are designated as 23 and 24, respectively
            ## So we need to check if primers are designed for human genome and chromosomes may be 23 and 24 instead of X and Y
            if species=='human':
                if chrom=='23': chromName='X'
                elif chrom=='24': chromName='Y'
                else: chromName=str(chrom)
            else:
                chromName=str(chrom)
            lines=pysam.faidx(wgref,'chr'+chromName+':'+str(chrTargetSeqStart)+'-'+str(amplBlockEnd+(args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen)+1)).split('\n')
            if lines[1]=='':
                print('ERROR! Extracted sequence has no length')
                print(chromName,start-1-args.maxAmplLen,end+args.maxAmplLen)
                exit(1)
            regionSeq=''.join(lines[1:-1]).upper()
            targetRegion=str(amplBlockStart-args.minPrimerShift-chrTargetSeqStart)+','+str(amplBlockEnd-amplBlockStart+1+2*args.minPrimerShift)
            targetRegionStart=amplBlockStart-args.minPrimerShift-chrTargetSeqStart
            targetRegionEnd=targetRegionStart+(amplBlockEnd-amplBlockStart+1+2*args.minPrimerShift)-1
            amplToStartCoord[amplName]=chrTargetSeqStart
            seqTags={}
            seqTags['SEQUENCE_TEMPLATE']=str(regionSeq)
            seqTags['SEQUENCE_ID']='NGS_primerplex_'+amplName
            # We create 3 sets of primer pairs, that will be located in three variants:
            ## --------------------Mut--------------------
            ## 1) -=======================-------------------
            ## 2) -------------------=======================-
            ## 3) Without forcing primer location
            for i in range(3):
                if i==0:
                    primerPairOkRegion=[[0,targetRegionStart-1,targetRegionEnd+1,len(regionSeq)-targetRegionEnd-1]]
                elif i==1:
                    # Force location of the right primer to the close proximity to target
                    primerPairOkRegion=[[0,targetRegionStart-1,targetRegionEnd+1,args.maxPrimerLen+2]]
                else:
                    # Force location of the left primer to the close proximity to target
                    primerPairOkRegion=[[targetRegionStart-args.maxPrimerLen-2,args.maxPrimerLen+2,targetRegionEnd+1,len(regionSeq)-targetRegionEnd-1]]
                seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                primer3Params[amplName].append([seqTags.copy(),primerTags.copy()])
        return(primer3Params,amplToChrom,amplToStartCoord)
    else:
        primerTags['PRIMER_PRODUCT_SIZE_RANGE']=[[args.minAmplLen,args.maxAmplLen]]
        primerTags['PRIMER_PRODUCT_OPT_SIZE']=args.optAmplLen
        # prevEnd stores coordinate of the previous position
        ## If currently processed position is more than previous end by one position,
        ## we need to make only primers that are near to the current position
        prevEnd=0
        for chrom,regions in pointRegions.items():
            for region in sorted(regions.values(),key=itemgetter(1)):
                chrom,start,end,curRegionName=region
                if regionNameToPrimerType is None:
                    primerTags['PRIMER_PICK_LEFT_PRIMER']=1
                    primerTags['PRIMER_PICK_RIGHT_PRIMER']=1
                elif curRegionName in regionNameToPrimerType.keys():
                    if regionNameToPrimerType[curRegionName]=='L':
                        primerTags['PRIMER_PICK_LEFT_PRIMER']=1
                        primerTags['PRIMER_PICK_RIGHT_PRIMER']=0
                    elif regionNameToPrimerType[curRegionName]=='R':
                        primerTags['PRIMER_PICK_LEFT_PRIMER']=0
                        primerTags['PRIMER_PICK_RIGHT_PRIMER']=1
                    elif regionNameToPrimerType[curRegionName]=='B':
                        primerTags['PRIMER_PICK_LEFT_PRIMER']=1
                        primerTags['PRIMER_PICK_RIGHT_PRIMER']=1
                    else:
                        print('ERROR! Unknown type of primers is necessary to be designed for the following region:')
                        print(curRegionName)
                        print('This value can be only "0" (only left primer), "1" (only right primer), "2" (both primers) or nothing (both primers)')
                        exit(1)
                # We do not need to split regions onto several blocks
                ## because previously we splited it onto point positions
                primer3Params[curRegionName]=[]
                # In COSMIC database chromosome X and Y are designated as 23 and 24, respectively
                ## So we need to check if primers are designed for human genome and chromosomes may be 23 and 24 instead of X and Y
                if species=='human':
                    if chrom=='23': chromName='X'
                    elif chrom=='24': chromName='Y'
                    else: chromName=str(chrom)
                else:
                    chromName=str(chrom)
                lines=pysam.faidx(wgref,'chr'+chromName+':'+str(start-args.maxAmplLen)+'-'+str(end+args.maxAmplLen)).split('\n')
                if lines[1]=='':
                    print('ERROR! Extracted sequence has no length')
                    print(chromName,start-1-args.maxAmplLen,end+args.maxAmplLen)
                    exit(1)
                regionSeq=''.join(lines[1:-1]).upper()
                seqTags={}
                seqTags['SEQUENCE_TEMPLATE']=str(regionSeq)
                seqTags['SEQUENCE_ID']='NGS_primerplex_'+curRegionName
                if start==prevEnd+1:                    
                    primerPairOkRegion=[[args.maxAmplLen-args.maxPrimerLen-2,args.maxPrimerLen+2,args.maxAmplLen+1+end-start,len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                    seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                    primer3Params[curRegionName].append([seqTags.copy(),primerTags.copy()])
                else:
                    # We create 3 sets of primer pairs, that will be located in three variants:
                    ## --------------------Mut--------------------
                    ## 1) -=======================-------------------
                    ## 2) -------------------=======================-
                    ## 3) Without forcing primer location
                    for i in range(3):
                        # Extracting reference sequence
                        if i==0:
                            # Force location of the right primer to the close proximity to target
                            primerPairOkRegion=[[0,args.maxAmplLen-1,args.maxAmplLen+1+end-start,args.maxPrimerLen+2]]
                        elif i==1:
                            # Force location of the left primer to the close proximity to target
                            primerPairOkRegion=[[args.maxAmplLen-args.maxPrimerLen-2,args.maxPrimerLen+2,args.maxAmplLen+1+end-start,len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                        else:
                            primerPairOkRegion=[[0,args.maxAmplLen-1,args.maxAmplLen+1+end-start,len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                        seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                        primer3Params[curRegionName].append([seqTags.copy(),primerTags.copy()])
                prevEnd=end
        return(primer3Params)

def constructInternalPrimers(primer3Params,regionNameToChrom,args,regionsCoords=None,allRegions=None,primersInfo=None,primersInfoByChrom=None,amplNames=None,primersToAmplNames=None):
    # chrom is string
    p=ThreadPool(args.threads)
    # Dictionary for storing primers' info
    if primersInfo==None:
        primersInfo={}
    # Dictionary for storing primers' info but primers are splitted by chromosome location
    if primersInfoByChrom==None:
        primersInfoByChrom={}
    # Dictionary for making unique amplicon names
    if amplNames==None:
        amplNames={}
    # Dictionary for converting primers' pairs to amplicon names
    if primersToAmplNames==None:
        primersToAmplNames={}
    # Constructing primers for each region with primer3
    totalPrimersNum=0
    totalDifPrimersNum=0
    regionsWithoutPrimers=[]
    results=[]
    wholeWork=len(primer3Params.keys())
    for i,(regionName,inputParams) in enumerate(primer3Params.items()):
        for inputParam in inputParams:
            results.append(p.apply_async(runPrimer3,(regionName,inputParam,False,args.autoAdjust)))
        showPercWork(i+1,wholeWork)
    print()
    doneWork=0
    wholeWork=len(results)
    primerDesignExplains={}
    for res in results:
        doneWork+=1
        showPercWork(doneWork,wholeWork)
        curRegionName,primerSeqs,primersCoords,primerTms,amplLens,amplScores,primer3File,designOutput=res.get()
        if curRegionName not in amplNames.keys():
            amplNames[curRegionName]=[]
        if primerSeqs==None:
            explanation=''
            if 'PRIMER_PAIR_EXPLAIN' in designOutput.keys():
                explanation=designOutput['PRIMER_PAIR_EXPLAIN']
            elif 'PRIMER_RIGHT_EXPLAIN' in designOutput.keys():
                explanation=designOutput['PRIMER_RIGHT_EXPLAIN']
            elif 'PRIMER_LEFT_EXPLAIN' in designOutput.keys():
                explanation=designOutput['PRIMER_LEFT_EXPLAIN']
            if curRegionName not in primerDesignExplains:
                primerDesignExplains[curRegionName]=[explanation]
            else:
                primerDesignExplains[curRegionName].append(explanation)
            continue
        chrom=regionNameToChrom[curRegionName]
        if int(chrom) not in primersInfoByChrom.keys():
            primersInfoByChrom[int(chrom)]={}
        # Extract start and end of target region
        try:
            start,end=allRegions[chrom][curRegionName][1:3]
        except KeyError:
            print('ERROR!',allRegions.keys())
            exit(1)
        # Go through each pair of primers and save them and info about them
        for i in range(int(len(primerSeqs)/2)):
            totalPrimersNum+=i+1
            # If this pair of primers has been already processed, then we save all possible info about it
            ## Including information that this primers cover other input positions
            if '_'.join(primerSeqs[2*i:2*i+2]) in primersInfo.keys():
                continue
            # Calculating coordinates of primers on a chromosome
            try:
                primersCoords[2*i][0]=primersCoords[2*i][0]+start-args.maxAmplLen # We do not substract 1, because we want to get real coordinate, not number of symbol in coordinate
            except TypeError:
                print('ERROR:',primersCoords)
                exit(1)
            # We substract args.maxAmplLen because it is an addiotional part of chromosome that we wrote to primer3 input file as target sequence
            # primersCoords[2*i+1][0] is a right coordinate of primer
            ## If we constructed only left primer
            if primersCoords[2*i+1][0]==0:
                # We set coordinate of the right primer as the most right position of this amplicon
                primersCoords[2*i+1][0]=start+args.maxAmplLen
            else:
                primersCoords[2*i+1][0]=primersCoords[2*i+1][0]+start-args.maxAmplLen
            # Calculate coordinates of amplified block that exclude primers
            amplBlockStart=primersCoords[2*i][0]+primersCoords[2*i][1]
            amplBlockEnd=primersCoords[2*i+1][0]-primersCoords[2*i+1][1]
            # Save info about this pair
            primersInfo['_'.join(primerSeqs[2*i:2*i+2])]=[primersCoords[2*i:2*i+2],primerTms[2*i:2*i+2],amplLens[i],amplScores[i],chrom]
            if '_'.join(primerSeqs[2*i:2*i+2]) not in primersInfoByChrom[int(chrom)].keys():
                primersInfoByChrom[int(chrom)]['_'.join(primerSeqs[2*i:2*i+2])]=primersCoords[2*i:2*i+2]
            # Check if this pair of primers covers also other target input regions
            coveredPointRegions=checkThisPrimerPairForCoveringOtherInputRegions(allRegions[chrom],amplBlockStart,amplBlockEnd)
            for curRegionName in coveredPointRegions:
                # Create name for this pair of primers and save it
                if curRegionName not in amplNames.keys():
                    amplNames[curRegionName]=[curRegionName+'_1']
                elif len(amplNames[curRegionName])==0:
                    amplNames[curRegionName].append(curRegionName+'_1')
                else:
                    prevNum=amplNames[curRegionName][-1].split('_')[-1]
                    amplNames[curRegionName].append(curRegionName+'_'+str(int(prevNum)+1))
                curAmplName=amplNames[curRegionName][-1]
                # Each pair of primers can cover several input positions
                if '_'.join(primerSeqs[2*i:2*i+2]) not in primersToAmplNames.keys():
                    primersToAmplNames['_'.join(primerSeqs[2*i:2*i+2])]=[curAmplName]
                else:
                    primersToAmplNames['_'.join(primerSeqs[2*i:2*i+2])].append(curAmplName)
            totalDifPrimersNum+=1
            ## TO DO:
            ### Check if this pair of primers has occured again due to genome repeat
    print()
    p.close()
    p.join()
    writeDraftPrimers(primersInfo,args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers.xls')
    regionsWithoutPrimers=[]
    for curRegionName,ampls in amplNames.items():
        if len(ampls)==0:
            regionsWithoutPrimers.append(curRegionName)
    if len(regionsWithoutPrimers)>0:
        print(' # WARNING! For ',len(regionsWithoutPrimers),'regions primers could not be designed with the defined parameters. Here are these regions:')
        logger.warn(' # WARNING! For '+str(len(regionsWithoutPrimers))+' regions primers could not be designed with the defined parameters. Here are these regions:')
        for regionWithoutPrimer in regionsWithoutPrimers:
            if args.skipUndesigned:
                regionsCoords[int(regionNameToChrom[regionWithoutPrimer])].remove(allRegions[regionNameToChrom[regionWithoutPrimer]][regionWithoutPrimer][1])
                if len(regionsCoords[int(regionNameToChrom[regionWithoutPrimer])])==0:
                    regionsCoords.pop(int(regionNameToChrom[regionWithoutPrimer]))
                allRegions[regionNameToChrom[regionWithoutPrimer]].pop(regionWithoutPrimer)
                regionNameToChrom.pop(regionWithoutPrimer)                
            print('   '+regionWithoutPrimer)
            logger.info('   '+regionWithoutPrimer)
            print(primer3Params[regionWithoutPrimer][0][0]['SEQUENCE_TEMPLATE'])
            logger.info(primer3Params[regionWithoutPrimer][0][0]['SEQUENCE_TEMPLATE'])
            for explanation in primerDesignExplains[regionWithoutPrimer]:
                print(explanation)
                logger.info(explanation)
        if not args.skipUndesigned:
            print(' You should use less stringent parameters')
            logger.info(' You should use less stringent parameters')
            exit(1)
    print(' # Total number of constructed primers:',totalPrimersNum)
    print(' # Total number of different constructed primers:',totalDifPrimersNum)
    logger.info(' # Total number of constructed primers: '+str(totalPrimersNum))
    logger.info(' # Total number of different constructed primers: '+str(totalDifPrimersNum))
    return(primersInfo,primersInfoByChrom,amplNames,primersToAmplNames,regionsCoords,regionNameToChrom,allRegions)    

def runPrimer3(regionName,inputParams,extPrimer,autoAdjust=False):
    extPrimerTryNum=0
    seqTags,primerTags=inputParams
    # Variable that stores, if we tried to set less stringent parameters for:
    ## polyN, GC-content of primers, length of primers,Tms of primers
    if autoAdjust:
        triedSofterParameters=[False,False,False,False,False]
    else:
        triedSofterParameters=[True,True,True,True,True]
    softerParametersNames=['polyN',"primers 3'-end GC-content",'primers GC-content','length of primers','Tm of primers']
    pParams=['PRIMER_MAX_POLY_X','PRIMER_MAX_END_GC','PRIMER_MIN_GC','PRIMER_MIN_SIZE','PRIMER_MIN_TM']
    while(True):
        try:
            out=primer3.designPrimers(seqTags,primerTags)
        except OSError as e:
            print('ERROR! '+str(e))
            print(seqTags,primerTags)
            logger.error('ERROR! '+str(e))
            logger.error(seqTags,primerTags)
            exit(1)
        numReturned=out['PRIMER_PAIR_NUM_RETURNED']
        if (primerTags['PRIMER_PICK_LEFT_PRIMER']==1 and
            primerTags['PRIMER_PICK_RIGHT_PRIMER']==1 and
            numReturned>0):
            # Check that all found primers has necessary GC-content of the 3'-end
            ## And all primers do not form hairpin with 5'-overhang
            minEndGC=primerTags['PRIMER_MIN_END_GC']
            primers=[]
            primersPoses=[]
            primersTms=[]
            primersProductSizes=[]
            primersProductPenalty=[]
            for i in range(numReturned):
                primers.append(out['PRIMER_LEFT_'+str(i)+'_SEQUENCE'].upper())
                primers.append(out['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'].upper())
                primersPoses.append(list(out['PRIMER_LEFT_'+str(i)]))
                primersPoses.append(list(out['PRIMER_RIGHT_'+str(i)]))
                primersTms.append(round(out['PRIMER_LEFT_'+str(i)+'_TM'],0))
                primersTms.append(round(out['PRIMER_RIGHT_'+str(i)+'_TM'],0))
                primersProductSizes.append(out['PRIMER_PAIR_'+str(i)+'_PRODUCT_SIZE'])
                primersProductPenalty.append(out['PRIMER_PAIR_'+str(i)+'_PENALTY'])
            # If after filtering no primers left and setting less stringent parameters did not solve the problem
            ## In the second iteration we do not filter primers by minimal GC-content of 3'-ends
##            if False in triedSofterParameters:
            i=0
            leftEnd3_rc=str(Seq.Seq(primers[2*i][-4:]).reverse_complement())
            rightEnd3_rc=str(Seq.Seq(primers[2*i+1][-4:]).reverse_complement())
            leftHairpin=primer3.calcHairpin(primers[2*i]).dg/1000
            rightHairpin=primer3.calcHairpin(primers[2*i+1]).dg/1000
            while(i<int(len(primers)/2) and len(primers)>0):
                if (primers[2*i][-5:].count('G')+primers[2*i][-5:].count('C')<minEndGC
                    or primers[2*i+1][-5:].count('G')+primers[2*i+1][-5:].count('C')<minEndGC
                    or (leftEnd3_rc in primers[2*i][1:-4-3] and leftHairpin<-2)
                    or (rightEnd3_rc in primers[2*i+1][1:-4-3] and rightHairpin<-2)):
                    primers.pop(2*i); primers.pop(2*i)
                    primersPoses.pop(2*i); primersPoses.pop(2*i)
                    primersTms.pop(2*i); primersTms.pop(2*i)
                    primersProductSizes.pop(i)
                    primersProductPenalty.pop(i)
                    numReturned-=1
                else:
                    i+=1
        elif (primerTags['PRIMER_PICK_LEFT_PRIMER']==1 and
              out['PRIMER_LEFT_NUM_RETURNED']>0):
            numReturned=out['PRIMER_LEFT_NUM_RETURNED']
            minEndGC=primerTags['PRIMER_MIN_END_GC']
            primers=[]
            primersPoses=[]
            primersTms=[]
            primersProductSizes=[]
            primersProductPenalty=[]
            for i in range(numReturned):
                primers.append(out['PRIMER_LEFT_'+str(i)+'_SEQUENCE'].upper())
                primers.append('')
                primersPoses.append(list(out['PRIMER_LEFT_'+str(i)]))
                primersPoses.append([0,0])
                primersTms.append(round(out['PRIMER_LEFT_'+str(i)+'_TM'],0))
                primersTms.append(0)
                primersProductSizes.append(0)
                primersProductPenalty.append(out['PRIMER_LEFT_'+str(i)+'_PENALTY'])
            i=0
            leftEnd3_rc=str(Seq.Seq(primers[2*i][-4:]).reverse_complement())
            leftHairpin=primer3.calcHairpin(primers[2*i]).dg/1000
            while(i<int(len(primers)/2) and len(primers)>0):
                if (primers[2*i][-5:].count('G')+primers[2*i][-5:].count('C')<minEndGC
                    or (leftEnd3_rc in primers[2*i][1:-4-3] and leftHairpin<-2)):
                    primers.pop(2*i); primers.pop(2*i)
                    primersPoses.pop(2*i); primersPoses.pop(2*i)
                    primersTms.pop(2*i); primersTms.pop(2*i)
                    primersProductSizes.pop(i)
                    primersProductPenalty.pop(i)
                    numReturned-=1
                else:
                    i+=1
        elif (primerTags['PRIMER_PICK_RIGHT_PRIMER']==1 and
              out['PRIMER_RIGHT_NUM_RETURNED']>0):
            numReturned=out['PRIMER_RIGHT_NUM_RETURNED']
            minEndGC=primerTags['PRIMER_MIN_END_GC']
            primers=[]
            primersPoses=[]
            primersTms=[]
            primersProductSizes=[]
            primersProductPenalty=[]
            for i in range(numReturned):
                primers.append('')
                primers.append(out['PRIMER_RIGHT_'+str(i)+'_SEQUENCE'].upper())
                primersPoses.append([0,0])
                primersPoses.append(list(out['PRIMER_RIGHT_'+str(i)]))
                primersTms.append(0)
                primersTms.append(round(out['PRIMER_RIGHT_'+str(i)+'_TM'],0))
                primersProductSizes.append(0)
                primersProductPenalty.append(out['PRIMER_RIGHT_'+str(i)+'_PENALTY'])
            i=0
            rightEnd3_rc=str(Seq.Seq(primers[2*i+1][-4:]).reverse_complement())
            rightHairpin=primer3.calcHairpin(primers[2*i+1]).dg/1000
            while(i<int(len(primers)/2) and len(primers)>0):
                if (primers[2*i+1][-5:].count('G')+primers[2*i+1][-5:].count('C')<minEndGC
                    or (rightEnd3_rc in primers[2*i+1][1:-4-3] and rightHairpin<-2)):
                    primers.pop(2*i); primers.pop(2*i)
                    primersPoses.pop(2*i); primersPoses.pop(2*i)
                    primersTms.pop(2*i); primersTms.pop(2*i)
                    primersProductSizes.pop(i)
                    primersProductPenalty.pop(i)
                    numReturned-=1
                else:
                    i+=1
        if numReturned==0:
            # If there is some parameter that we haven't tried to set softer
            if False in triedSofterParameters:
                # Found this parameter
                paramNum=triedSofterParameters.index(False)
                if paramNum<=1:
                    curValue=primerTags[pParams[paramNum]]
                    newValue=min(5,curValue+1)
                    primerTags[pParams[paramNum]]=newValue
                    triedSofterParameters[paramNum]=True
                    continue
                else:
                    pParamMin=pParams[paramNum]
                    pParamMax=pParams[paramNum].replace('MIN','MAX')
                    curValueMin=primerTags[pParamMin]
                    curValueMax=primerTags[pParamMax]
                    if paramNum==2:
                        newValueMin=max(0,curValueMin-5)
                        newValueMax=min(100,curValueMax+5)
                    elif paramNum==3:
                        newValueMin=max(15,curValueMin-1)
                        newValueMax=min(36,curValueMax+1)
                    elif paramNum==4:
                        newValueMin=max(40,curValueMin-1)
                        newValueMax=min(95,curValueMax+1)
                    primerTags[pParamMin]=newValueMin
                    primerTags[pParamMax]=newValueMax
                    triedSofterParameters[paramNum]=True
                    continue
            elif seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'][0][0]!=0:
                seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'][0][0]=0
                continue
            else:
                return(regionName,None,None,None,None,None,inputParams,out)
        else:
            return(regionName,primers,primersPoses,primersTms,primersProductSizes,primersProductPenalty,inputParams,out)

def writeDraftPrimers(primersInfo,rFile,goodPrimers=None,external=False):
    # Good primers contains list of primers that correspond some defined parameters, e.g. specificity
    # worksheet #0 - internal primers
    # worksheet #1 - external primers
    oldRows=[]
    if os.path.exists(rFile):
        wb=xlrd.open_workbook(rFile)
        ws=None
        if external:
            ws=wb.sheet_by_index(0)
        elif wb.nsheets>1:
            ws=wb.sheet_by_index(1)
        if ws:
            for i in range(ws.nrows):
                if i==0:
                    continue
                row=ws.row_values(i)
                oldRows.append(row)
    wbw=xls.Workbook(rFile)
    if external:
        wsw1=wbw.add_worksheet('Draft_Internal_Primers')
        wsw2=wbw.add_worksheet('Draft_External_Primers')
    else:
        wsw2=wbw.add_worksheet('Draft_Internal_Primers')
        wsw1=wbw.add_worksheet('Draft_External_Primers')
    wsw1.write_row(0,0,['Primer_Pair','Left_Primer_Start','Left_Primer_Length',
                        'Right_Primer_End','Right_Primer_Length',
                        'Left_Primer_Tm','Right_Primer_Tm',
                        'Amplicon_Length','Primers_Score','Chrom'])
    rowNum=1
    for row in oldRows:
        wsw1.write_row(rowNum,0,row)
        rowNum+=1
    wsw2.write_row(0,0,['Primer_Pair','Left_Primer_Start','Left_Primer_Length',
                       'Right_Primer_End','Right_Primer_Length',
                       'Left_Primer_Tm','Right_Primer_Tm',
                       'Amplicon_Length','Primers_Score','Chrom'])
    rowNum=1
    for primerPair,info in primersInfo.items():
        if (goodPrimers==None or
            goodPrimers and primerPair in goodPrimers):
            wsw2.write_row(rowNum,0,[primerPair,info[0][0][0],info[0][0][1],
                                    info[0][1][0],info[0][1][1],
                                    info[1][0],info[1][1],
                                    info[2],info[3],info[4]])
            rowNum+=1
    wbw.close()

def readDraftPrimers(draftFile,external=False):
    wb=xlrd.open_workbook(draftFile)
    if external:
        ws=wb.sheet_by_index(1)
    else:
        ws=wb.sheet_by_index(0)
    # Read info about designed primers
    primersInfo={}
    primersInfoByChrom={}
    for i in range(ws.nrows):
        if i==0:
            continue
        row=ws.row_values(i)
        ##    [primersCoords[2*i:2*i+2],primerTms[2*i:2*i+2],amplLens[i],amplScores[i],chrom]
        primersInfo[row[0]]=[[[int(row[1]),int(row[2])],[int(row[3]),int(row[4])]],
                             [int(row[5]),int(row[6])],
                             int(row[7]),int(row[8]),str(row[9])]
        if getChrNum(row[9]) not in primersInfoByChrom.keys():
            primersInfoByChrom[getChrNum(row[9])]={row[0]:[[int(row[1]),int(row[2])],[int(row[3]),int(row[4])]]}
        elif row[0] not in primersInfoByChrom[getChrNum(row[9])].keys():
            primersInfoByChrom[getChrNum(row[9])][row[0]]=[[int(row[1]),int(row[2])],[int(row[3]),int(row[4])]]
    return(primersInfo,primersInfoByChrom)

def getChrNum(chrom):
    if 'chr' in str(chrom):
        chrom=chrom.replace('chr','')
    if chrom=='X' or chrom=='x':
        chrom=23
    elif chrom=='Y' or chrom=='y':
        chrom=24
    else:
        chrom=int(chrom)
    return(chrom)

def getRegionsUncoveredByDraftPrimers(allRegions,primersInfoByChrom):
    # allRegions[chrom]={curRegionName:[chrom,i,i+endShift,curRegionName]}
    # primersInfoByChrom[int(chrom)]['_'.join(primerSeqs[2*i:2*i+2])]=primersCoords[2*i:2*i+2]
    amplNames={}
    primersToAmplNames={}
    uncoveredRegions=allRegions.copy()
    for chrom,coords in primersInfoByChrom.items():
        coordToRegionName={}
        for regionCoords in allRegions[str(chrom)].values():
            coordToRegionName[regionCoords[1]]=regionCoords[3]
        for primerPair,coord in coords.items():
            for i in range(coord[0][0]+coord[0][1],coord[1][0]-coord[1][1]+1):
                if i in coordToRegionName.keys():
                    curRegionName=coordToRegionName[i]
                    # Create name for this pair of primers and save it
                    if curRegionName not in amplNames.keys():
                        amplNames[curRegionName]=[curRegionName+'_1']
                    else:
                        prevNum=amplNames[curRegionName][-1].split('_')[-1]
                        amplNames[curRegionName].append(curRegionName+'_'+str(int(prevNum)+1))
                    curAmplName=amplNames[curRegionName][-1]
                    # Each pair of primers can cover several input positions
                    if primerPair not in primersToAmplNames.keys():
                        primersToAmplNames[primerPair]=[curAmplName]
                    else:
                        primersToAmplNames[primerPair].append(curAmplName)
                    if (curRegionName in uncoveredRegions[str(chrom)].keys() and
                        allRegions[str(chrom)][curRegionName][1]==allRegions[str(chrom)][curRegionName][2]):
                        uncoveredRegions[str(chrom)].pop(curRegionName)
                    elif (curRegionName in uncoveredRegions[str(chrom)].keys() and
                          coord[0][0]+coord[0][1]<=allRegions[str(chrom)][curRegionName][1] and
                          coord[1][0]-coord[1][1]>=allRegions[str(chrom)][curRegionName][2]):
                        uncoveredRegions[str(chrom)].pop(curRegionName)
                        break
        if len(uncoveredRegions[str(chrom)])==0:
            uncoveredRegions.pop(str(chrom))
    return(uncoveredRegions,amplNames,primersToAmplNames)

def getRegionsUncoveredByDraftExternalPrimers(primersInfo,primersInfoByChrom,outputInternalPrimers,args):
    # primersInfo[primer1_primer2]=[[[leftStart,primerLength],[rightEnd,rightLength]],
                                  # [leftTm,rightTm],
                                  # amplLen,score,chromosome]
    # primersInfoByChrom[int(chrom)]['_'.join(primerSeqs[2*i:2*i+2])]=primersCoords[2*i:2*i+2]
    # outputInternalPrimers[amplName]=[internalPrimers[0],internalPrimers[1],amplName,chrom,
                                     # amplStart,amplEnd,amplLen,amplBlockStart,
                                     # amplBlockEnd,leftPrimerTm,rightPrimerTm,len(internalPrimers[0]),
                                     # len(internalPrimers[1]),leftGC,rightGC]
    # outputExternalPrimers[curRegionName]=[[primerSeqs[2*k+0],primerSeqs[2*k+1],curRegionName+'_ext',chrom,start,end,
                                           # end-start+1,start+primersCoords[2*k][1],end-primersCoords[2*k+1][1],
                                           # primerTms[2*k+0],primerTms[2*k+1],len(primerSeqs[2*k+0]),len(primerSeqs[2*k+1]),
                                           # leftGC,rightGC,outputInternalPrimers[curRegionName][7]-(start+primersCoords[0][1])+1,
                                           # end-primersCoords[1][1]-outputInternalPrimers[curRegionName][8]+1,extendedAmplSeq]]
    outputExternalPrimers={}
    uncoveredInternalPrimers=outputInternalPrimers.copy()
    amplToChrom={}
    amplToStartCoord={}
    refFa=pysam.FastaFile(args.wholeGenomeRef)
    for chrom,coords in primersInfoByChrom.items():
        for amplName,parameters in sorted(outputInternalPrimers.items(),
                                          key=lambda item:item[1][3]):
            if str(chrom)!=parameters[3]:
                continue
            elif str(chrom)>parameters[3]:
                break
            for primerPair,coord in coords.items():
                # If current external primers overlap with internal no more than defined shift
                if (parameters[7]-1-(coord[0][0]+coord[0][1]-1)>=args.minPrimerShift and
                    coord[1][0]-coord[1][1]+1-parameters[8]+1>=args.minPrimerShift):
                    leftPrimer,rightPrimer=primerPair.split('_')
                    if len(leftPrimer)>0:
                        leftGC=round((leftPrimer.count('G')+leftPrimer.count('C'))*100/len(leftPrimer),2)
                    else:
                        leftGC=0
                    if len(rightPrimer)>0:
                        rightGC=round((rightPrimer.count('G')+rightPrimer.count('C'))*100/len(rightPrimer),2)
                    else:
                        rightGC=0
                    seq=extractGenomeSeq(refFa,chrom,coord[0][0]-100,coord[1][0]+100)
                    outputExternalPrimers[parameters[2]]=[[leftPrimer,rightPrimer,parameters[2]+'_ext',
                                                          str(chrom),coord[0][0],coord[1][0],coord[1][0]-coord[0][0]+1,
                                                          coord[0][0]+coord[1][0],coord[1][0]-coord[1][1],
                                                          primersInfo[primerPair][1][0],primersInfo[primerPair][1][1],
                                                          len(leftPrimer),len(rightPrimer),leftGC,rightGC,
                                                          parameters[7]-1-(coord[0][0]+coord[0][1]-1),
                                                          coord[1][0]-coord[1][1]+1-parameters[8]+1,seq]]
                    uncoveredInternalPrimers.pop(amplName)
                    chrTargetSeqStart=parameters[9]-(args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen)
                    amplToChrom[parameters[2]]=str(chrom)
                    amplToStartCoord[parameters[2]]=chrTargetSeqStart
    return(outputExternalPrimers,uncoveredInternalPrimers,amplToChrom,amplToStartCoord)

def checkThisPrimerPairForCoveringOtherInputRegions(chromPointRegions,amplBlockStart,amplBlockEnd):
    coveredRegions=[]
    allStarts=[region[1] for region in chromPointRegions.values()]
    allStarts.append(amplBlockStart)
    allStarts.sort()
    for region in sorted(chromPointRegions.values(),key=itemgetter(1))[allStarts.index(amplBlockStart):]:
        chrom,start,end,curRegionName=region
        if start>amplBlockEnd:
            break
        elif start>=amplBlockStart and end<=amplBlockEnd:
            coveredRegions.append(curRegionName)
    return(coveredRegions)

def checkPrimersSpecificity(inputFileBase,primersInfo,wholeGenomeRef,runName,substNum=1,threads=2,maxNonSpecLen=100,maxPrimerNonspec=1000,external=False,varNum=''):
    # Dictionary for storing info about primers specificity by primers
    ## for checking specificity within one amplicon    
    primersNonSpecRegions={}
    # Dictionary for storing info about primers specificity by chromosome
    ## for checking specificity within one multiplex
    primersNonSpecRegionsByChrs={}
    # Creating fasta-file with all primers' sequences
    if external:
        seqFile=open(inputFileBase+'_NGS_primerplex'+runName+'_all_external_primers'+varNum+'_sequences.fa','w')
    else:
        seqFile=open(inputFileBase+'_NGS_primerplex'+runName+'_all_primers_sequences.fa','w')
    for key in primersInfo.keys():
        primers=key.split('_')
        for j,primer in enumerate(primers):
            if primer!='' and primer not in primersNonSpecRegions.keys():
                seqFile.write('\n'.join(['>'+primer,primer])+'\n')
                primersNonSpecRegions[primer]=None   
            elif primer not in primersNonSpecRegions.keys():
                primersNonSpecRegions[primer]=None                
    seqFile.close()
    if external:
        bwaResultFileName=inputFileBase+'_NGS_primerplex'+runName+'_all_external_primers'+varNum+'_sequences.bwa'
    else:
        bwaResultFileName=inputFileBase+'_NGS_primerplex'+runName+'_all_primers_sequences.bwa'
    print(' Running BWA...')
    logger.info(' Running BWA...')
    out=sp.check_output('bwa aln -N -n '+str(args.substNum)+' -t '+str(threads)+' '+wholeGenomeRef+' '+seqFile.name+' > '+bwaResultFileName+'.sai',shell=True,stderr=sp.STDOUT).decode('utf-8')
    out=sp.check_output('bwa samse -n 10000000 '+wholeGenomeRef+' '+bwaResultFileName+'.sai'+' '+seqFile.name+' > '+bwaResultFileName+'.sam',shell=True,stderr=sp.STDOUT).decode('utf-8')
    # Reading BWA output file
    samFile=pysam.AlignmentFile(bwaResultFileName+'.sam')
    # Process SAM-file strings in several threads
    p=ThreadPool(threads)
    print(' Processing SAM-file...')
    logger.info(' Processing SAM-file...')
    totalNumberNonspecificRegions=0
    results=[]
    refFa=pysam.FastaFile(wholeGenomeRef)
    for read in samFile.fetch():
        results.append(p.apply_async(readBwaFile,(read,maxPrimerNonspec,refFa,primersInfo)))
    wholeWork=len(results)
    doneWork=0
    unspecificPrimers={}
    for res in results:
        # res[0] is a primer name
        res=res.get()
        doneWork+=1
        showPercWork(doneWork,wholeWork)
        if res[1]==None or len(res[1])==0:
            primersNonSpecRegions[res[0]]=None
        else:
            primersNonSpecRegions[res[0]]={}
            totalNumberNonspecificRegions+=len(res[1])
            if len(res[1])>maxPrimerNonspec:
                unspecificPrimers[res[0]]=len(res[1])
                continue
            # res[1] has format: chrom, pos, length
            for region in res[1]:
                if region[0] not in primersNonSpecRegions[res[0]].keys():
                    primersNonSpecRegions[res[0]][region[0]]=[region[1:]]
                else:
                    primersNonSpecRegions[res[0]][region[0]].append(region[1:])
                if region[0] not in primersNonSpecRegionsByChrs.keys():
                    primersNonSpecRegionsByChrs[region[0]]={res[0]:[region[1:]]}
                elif res[0] not in primersNonSpecRegionsByChrs[region[0]].keys():
                    primersNonSpecRegionsByChrs[region[0]][res[0]]=[region[1:]]
                else:
                    primersNonSpecRegionsByChrs[region[0]][res[0]].append(region[1:])
    p.close()
    p.join()
    refFa.close()
    print('\n # Total number of nonspecific regions:',totalNumberNonspecificRegions)
    logger.info(' # Total number of nonspecific regions: '+str(totalNumberNonspecificRegions))
    # Now we go through all primer pairs and check them for nonspecific amplicons within one pair of primers
    print(' Searching for nonspecific amplicons that are formed by designed primer pairs...')
    logger.info(' Searching for nonspecific amplicons that are formed by designed primer pairs...')
    # File for storing statistics of nonspecific regions about each primer
    if varNum!='':
        varNum='_'+varNum
    if not external:
        # This file includes specificity information about all primers including those one that were not included to any output combination
        wbw_spec=xls.Workbook(inputFileBase+'_NGS_primerplex'+runName+'_all_internal_primers_specificity.xls')
    else:
        wbw_spec=xls.Workbook(inputFileBase+'_NGS_primerplex'+runName+'_external_primers'+varNum+'_specificity.xls')
    wsw_spec=wbw_spec.add_worksheet('Primers_specificity')
    wsw_spec.write_row(0,0,['Primer','Number_of_Nonspecific_regions'])
    colsWidth=[30,30]
    for k,colsWidth in enumerate(colsWidth):
        wsw_spec.set_column(k,k,colsWidth)
    specificPrimers=[]
    rowNum=1
    for primerPair in primersInfo.keys():
        primerName1,primerName2=primerPair.split('_')
        hasNonspecAmpl=False
        # If both primers have too many unspecific regions for annealing
        if primerName1 in unspecificPrimers.keys() and primerName2 in unspecificPrimers.keys():
            hasNonspecAmpl=True
            totalNonspecRegionsForPrimer1=unspecificPrimers[primerName1]
            totalNonspecRegionsForPrimer2=unspecificPrimers[primerName2]
        elif primerName1 in unspecificPrimers.keys():
            hasNonspecAmpl=True
            totalNonspecRegionsForPrimer1=unspecificPrimers[primerName1]
            totalNonspecRegionsForPrimer2='N/A'
        elif primerName2 in unspecificPrimers.keys():
            hasNonspecAmpl=True
            totalNonspecRegionsForPrimer2=unspecificPrimers[primerName2]
            totalNonspecRegionsForPrimer1='N/A'
        else:
            totalNonspecRegionsForPrimer1=0
            totalNonspecRegionsForPrimer2=0
            regions1=primersNonSpecRegions[primerName1]
            regions2=primersNonSpecRegions[primerName1]
            if regions1!=None and regions2==None:
                for chrRegions in regions1.values():
                    totalNonspecRegionsForPrimer1+=len(chrRegions)
            elif regions2!=None and regions1==None:
                for chrRegions in regions2.values():
                    totalNonspecRegionsForPrimer2+=len(chrRegions)
            elif regions1!=None and regions2!=None:
                # We go through only one of dict because we are going to search nonspecific amplicons only 
                ## on chromosomes that are presented in both dicts
                for chrom,chrRegions in regions1.items():
                    totalNonspecRegionsForPrimer1+=len(chrRegions)
                    if chrom not in regions2.keys(): continue
                    for reg1 in chrRegions:
                        for reg2 in regions2[chrom]:
                            # Strands should be different for two nonspecific regions of primers
                            if ((reg1[0]>0 and reg2[0]<0 and 0<reg2[1]+reg2[2]-reg1[1]<=maxNonSpecLen) or
                                (reg1[0]<0 and reg2[0]>0 and 0<reg1[1]+reg1[2]-reg2[1]<=maxNonSpecLen)):
                                # This means that this pair of primers form nonspecific amplicons
                                ## So we should remove it
                                hasNonspecAmpl=True
                                break
                        if hasNonspecAmpl: break
                    if hasNonspecAmpl: break
                for chrRegions in regions2.values():
                    totalNonspecRegionsForPrimer2+=len(chrRegions)
        if (not hasNonspecAmpl and (totalNonspecRegionsForPrimer1>maxPrimerNonspec or
                                   totalNonspecRegionsForPrimer2>maxPrimerNonspec)):
            hasNonspecAmpl=True
        elif not hasNonspecAmpl:
            specificPrimers.append(primerPair)
        wsw_spec.write_row(rowNum,0,[primerName1,totalNonspecRegionsForPrimer1])
        wsw_spec.write_row(rowNum+1,0,[primerName2,totalNonspecRegionsForPrimer2])
        rowNum+=2 
    wbw_spec.close()
    # each element of specificPrimers has a format: primerSeq1_primerSeq2
    return(specificPrimers,primersNonSpecRegionsByChrs)

def readBwaFile(read,maxPrimerNonspec,refFa,primersInfo=None):
    indelsPat=re.compile('(\d+)([ID])')
    matchPat=re.compile('(\d+)M')
    # Extract strand of the main match in genome
    if read.flag==0:
        strand=1
    elif read.flag==16:
        strand=-1
    else:
        print('ERROR! Unknown value of FLAG:',read.flag)
        print(read)
        exit(1)
    if primersInfo:
        primersName=read.qname#[:read.qname.rfind('_')]
        for primerPair in primersInfo.keys():
            if primersName in primerPair.split('_'):
                try:
                    targetRegion=['chr'+primersInfo[primerPair][4],int(primersInfo[primerPair][0][primerPair.split('_').index(primersName)][0])]
                except TypeError:
                    print('ERROR!',primersInfo[primerPair][4],primersInfo[primerPair][0],primerPair.split('_').index(primersName))
                    exit(1)
                break
    else:
        targetRegion=[]
    if [read.reference_name,strand*(read.pos+1)]!=targetRegion:
        nonSpecRegions=[','.join([read.reference_name,str(strand*(read.pos+1)),
                                  read.cigarstring,str(read.get_tag('NM'))])]
    else:
        nonSpecRegions=[]
    if read.has_tag('XA'):
        nonSpecRegions.extend(read.get_tag('XA').split(';')[:-1])
    qname='_'.join(read.qname.split('_')[:2]+read.qname.split('_')[-1:])
    if len(nonSpecRegions)>maxPrimerNonspec*2:
        return(qname,nonSpecRegions)
    primerNonSpecRegions=[]
    for region in nonSpecRegions: # read XA tag ends with ; so the last element is empty
        # We go through all these regions and check that 3'-nucleotide matches primer's 3'-end
        chrom,pos,cigar,subst=region.split(',')
        # Determine length of sequence of interest by parsing CIGAR
        ## The length depends only on the deletions but not mismatches nor insertions into reference genome
        ### So we count number of deletions
        indelsMatch=indelsPat.findall(cigar)
        matchMatch=matchPat.findall(cigar)
        if len(indelsMatch)==0:
            regionLen=int(matchMatch[0])
        else:
            sumDeletions=0
            sumMatch=0
            for m in indelsMatch:
                if m[1]=='D':
                    sumDeletions+=int(m[0])
            for m in matchMatch:
                sumMatch+=int(m)
            regionLen=sumMatch+sumDeletions
        attempts=0
        seq=None
        while(seq is None):
            try:
                seq=refFa.fetch(region=chrom+':'+str(abs(int(pos)))+'-'+str(abs(int(pos))+regionLen-1))
            except:
                seq=None
                attempts+=1
                if attempts>=10:                    
                    print('ERROR!',e)
                    logger.error(str(e))
                    print(refFa.filename)
                    logger.error(refFa.filename)
                    exit(1)
        # Remove 
        # Determine, which sequence we should take: forward or reverse-complement
        ## If primer is on + strand and found region is on opposite
        ### or primer is on - strand and found region is on the same
        #### we take reverse-complement
        if int(pos)>0:
            regStrand=1
        elif int(pos)<0:
            regStrand=-1
        if (strand>0 and int(pos)<0) or (strand<0 and int(pos)>0):
            seq=str(Seq.Seq(seq).reverse_complement())
        # If there is some insertions or deletion in found region
        if 'I' in cigar or 'D' in cigar:
            # We need to align its sequence with primer sequence
            align=pairwise2.align.globalxx(read.seq,seq)
            # Check that 3'-ends of primers are identical
            ## and one of two nucleotides before 3'-ends are identical, too
            if align[0][0][-1]==align[0][1][-1] and (align[0][0][-2]==align[0][1][-2] or
                                                     align[0][0][-3]==align[0][1][-3]):
                # Then we consider this region as a non-specific for this primer
                primerNonSpecRegions.append([chrom,regStrand,abs(int(pos)),regionLen])
        else:
            if seq[-1]==read.seq[-1] and (seq[-2]==read.seq[-2] or
                                         seq[-3]==read.seq[-3]):
                # Then we consider this region as a non-specific for this primer
                primerNonSpecRegions.append([chrom,regStrand,abs(int(pos)),regionLen])
    return(read.qname,primerNonSpecRegions)

def getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,maxNonSpecLen=100,threads=2):
    unspecificPrimers=set()
    p=ThreadPool(threads)
    results=[]
    for chrom,primers in primersNonSpecRegionsByChrs.items():
        sortedPrimers=sorted(primers.items())
        for i,(primer1,regions1) in enumerate(sortedPrimers[:-1]):
            results.append(p.apply_async(checkPrimerCoordinatesWithOtherPrimers,(primer1,regions1,sortedPrimers[i+1:],maxNonSpecLen)))
    allWork=len(results)
    for i,res in enumerate(results):
        res=res.get()
        unspecificPrimers.update(res)
        showPercWork(i+1,allWork)
    return(unspecificPrimers)

def checkPrimerCoordinatesWithOtherPrimers(primer1,regions1,primers2,maxNonSpecLen):
    unspecificPrimers=set()
    for reg1 in regions1:
        for primer2,regions2 in sorted(primers2):
            for reg2 in regions2:
                # Strands should be different for two nonspecific regions of primers
                if (reg1[0]>0 and reg2[0]<0 and 0<reg2[1]+reg2[2]-reg1[1]<=maxNonSpecLen or
                    reg1[0]<0 and reg2[0]>0 and 0<reg1[1]+reg1[2]-reg2[1]<=maxNonSpecLen):
                    # Extract only primer sequence that forms nonspecific product
                    fPrimer1=primer1
                    rPrimer2=primer2
                    unspecificPrimers.add('_'.join([fPrimer1,rPrimer2]))
    return(unspecificPrimers)

def removeBadPrimerPairs(primersInfoByChrom,goodPrimers,primersToAmplNames,amplNames):
    newPrimersInfoByChrom={}
    for chrom,primers in primersInfoByChrom.items():
        for primerPairName,primers in primers.items():
            if primerPairName in goodPrimers:
                if chrom not in newPrimersInfoByChrom.keys():
                    newPrimersInfoByChrom[chrom]={primerPairName:primers}
                elif primerPairName not in newPrimersInfoByChrom[chrom].keys():
                    newPrimersInfoByChrom[chrom][primerPairName]=primers
                else:
                    print('ERROR! Pair of primers is repeated in the primersInfoByChrom!')
                    print(chrom,primerPairName)
                    logger.error('Pair of primers is repeated in the primersInfoByChrom!')
                    logger.error(chrom,primerPairName)
                    exit(1)
            else:
                amplNamesToDelete=primersToAmplNames[primerPairName]
                for amplName in amplNamesToDelete:
                    curRegionName=amplName[:amplName.rfind('_')]
                    amplNames[curRegionName].remove(amplName)
    primersInfoByChrom=newPrimersInfoByChrom.copy()
    return(primersInfoByChrom,amplNames)

def checkThatAllInputRegionsCovered(amplNames,allRegions,regionNameToChrom,regionsCoords,filterMessage='by specificity',skipUndesigned=False):
    someInputRegionUncovered=False
    newAmplNames={}
    for curRegionName,ampls in amplNames.items():
        if len(ampls)==0:
            if skipUndesigned:
                print('WARNING! For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                logger.warn('For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                regionsCoords[int(regionNameToChrom[curRegionName])].remove(allRegions[regionNameToChrom[curRegionName]][curRegionName][1])
                if len(regionsCoords[int(regionNameToChrom[curRegionName])])==0:
                    regionsCoords.pop(int(regionNameToChrom[curRegionName]))
                allRegions[regionNameToChrom[curRegionName]].pop(curRegionName)
                regionNameToChrom.pop(curRegionName)
            else:
                newAmplNames[curRegionName]=ampls
                print('ERROR! For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                logger.error('For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                someInputRegionUncovered=True
        else:
            newAmplNames[curRegionName]=ampls
    if someInputRegionUncovered:
        exit(1)
    return(newAmplNames,allRegions,regionNameToChrom,regionsCoords)

def analyzePrimersForCrossingSNP(primersInfo,threads,localGenomeDB={},genomeVersion='hg19'):
    p=ThreadPool(threads)
    results=[]
    for primers,info in sorted(primersInfo.items(),
                               key=lambda item:(item[1][4],item[1][0][0][0])):
        primer1,primer2=primers.split('_')
        chrom=info[4]
        try:
            start1=info[0][0][0]; start2=info[0][1][0]-info[0][1][1]+1
            end1=info[0][0][0]+info[0][0][1]-1; end2=info[0][1][0]
        except IndexError as e:
            print('ERROR!',e)
            print(info)
            exit(1)
        strand1=1; strand2=-1
        if primer1!='':
            results.append(p.apply_async(checkPrimerForCrossingSNP,(primer1,chrom,start1,end1,strand1,localGenomeDB,args.snpFreq,genomeVersion)))
        if primer2!='':
            results.append(p.apply_async(checkPrimerForCrossingSNP,(primer2,chrom,start2,end2,strand2,localGenomeDB,args.snpFreq,genomeVersion)))
    primersCoveringSNPs=[]
    wholeWork=len(results)
    done=0
    for res in results:
        primer,result,localGenomeDB=res.get()
        if result:
            primersCoveringSNPs.append(primer)
        done+=1
        showPercWork(done,wholeWork)
    p.close()
    p.join()
    primerPairsNonCoveringSNPs=[]
    for primers,info in primersInfo.items():
        primer1,primer2=primers.split('_')
        if primer1 not in primersCoveringSNPs and primer2 not in primersCoveringSNPs:
            primerPairsNonCoveringSNPs.append(primers)
    print("\n # Number of primers that do not overlap with high-frequent SNPs by 3'-end: "+str(len(primerPairsNonCoveringSNPs)))
    logger.info("\n # Number of primers that do not overlap with high-frequent SNPs by 3'-end: "+str(len(primerPairsNonCoveringSNPs)))
    return(primerPairsNonCoveringSNPs,primersCoveringSNPs,localGenomeDB)    

# checkPrimerForCrossingSNP checks, if primer crosses some SNPs with high frequence in population
def checkPrimerForCrossingSNP(primer,chrom,start,end,strand,localGenomeDB={},freq=0.1,genomeVersion='hg19'):
    mv=myvariant.MyVariantInfo()
    # Previously, we only filtered out primers that crossed SNP by 3'-end
    # Now, we filter primers that cross SNP by any of its sequence
    for i,coord in enumerate(range(start,end+1)):
        if strand>0:
            nuc=primer[i]
##            coordStr='_'.join([str(chrom),str(coord)])
        else:
            nuc=primer[-(i+1)]
        coordStr='_'.join([str(chrom),str(coord)])
        if coordStr not in localGenomeDB.keys():
            try:
                if strand>0:
                    mvRes=mv.query('dbsnp.chrom:'+str(chrom)+' && dbsnp.'+genomeVersion+'.start:'+str(coord),fields='dbsnp')
                    coord=end
                else:
                    mvRes=mv.query('dbsnp.chrom:'+str(chrom)+' && dbsnp.'+genomeVersion+'.start:'+str(coord),fields='dbsnp')
                    coord=start
            except:
                print('ERROR! Could not check primer for crossing SNPs with the following parameters:')
                print(primer,chrom,start,end,strand)
                exit(1)
            hfSnpFound=False
            for hit in mvRes['hits']:
                if 'gmaf' not in hit['dbsnp'].keys(): continue
                alFreqs={}
                for al in hit['dbsnp']['alleles']:
                    if 'freq' not in al.keys(): continue
                    try:
                        alFreqs[al['allele']]=al['freq']
                    except KeyError:
                        print('ERROR!',mvRes); exit(1)
                ref=hit['dbsnp']['ref']
                alt=hit['dbsnp']['alt']
                if alt not in alFreqs.keys(): continue
                try:
                    if alFreqs[alt]>=freq:
                        hfSnpFound=True
                except KeyError:
                    print('ERROR of key:',alFreqs)
                    print(mvRes)
                    exit(1)
            if coordStr not in localGenomeDB.keys():
                localGenomeDB[coordStr]=hfSnpFound
            else:
                localGenomeDB[coordStr]=hfSnpFound
            if hfSnpFound:
                return(primer,hfSnpFound,localGenomeDB)
        else:
            return(primer,localGenomeDB[coordStr],localGenomeDB)
    return(primer,hfSnpFound,localGenomeDB)

# joinAmpliconsToAmplifiedBlocks joins neighbourhing or overlapping amplicons to blocks
# Minimal path is a path of one primer:
## [coord1,primer_name,coord2]
def joinAmpliconsToBlocks(chromRegionsCoords,chromPrimersInfoByChrom,maxAmplLen=100,chrom=None):
    # First, split all regions on this chromosome onto blocks, elements of which cannot be joined into one amplificated block
    blocks=[[chromRegionsCoords[0]]]
    coordToBlock={} # converts coordinate into the number of block 
    for coord in chromRegionsCoords[1:]:
        if coord>=blocks[-1][-1]+maxAmplLen:
            blocks.append([coord])
        else:
            blocks[-1].append(coord)
        coordToBlock[coord]=len(blocks)-1
    logger.info('Chromosome '+str(chrom)+' was splitted onto '+str(len(blocks))+' distinct blocks')
    firstNodes=[]
    lastNodes=[]
    blockGraphs=[]
    for block in blocks:
        # Making graphs for each chromosome
        blockGraphs.append(nx.Graph())
        # Get the first and the last input positions on this chromosome
        firstNodes.append(block[0])
        lastNodes.append(block[-1])
    # Add first mutation as the first node and search for primer pairs that cover it
    ## Also for each primer pair we search for primer pairs that cover the next input position
    ### For each edge we add weight, that describes farness from that edge
    for i,(primerPairName1,primers1) in enumerate(sorted(chromPrimersInfoByChrom.items(),key=lambda item:item[1][0][0])):
        primerName1,primerName2=primerPairName1.split('_')
        amplBlockStart1=primers1[0][0]+primers1[0][1]
        amplBlockEnd1=primers1[1][0]-primers1[1][1]
        # Get block number
        coords=chromRegionsCoords.copy()
        coords.append(amplBlockStart1)
        amplBlockStart1_num=sorted(coords).index(amplBlockStart1)
        if amplBlockStart1_num==0:
            blockNum=0
        elif amplBlockStart1_num==len(chromRegionsCoords):
            blockNum=len(blocks)-1
        else:
            nextMut=chromRegionsCoords[amplBlockStart1_num]
            blockNum=coordToBlock[nextMut]
        firstMut=False
        lastMut=False
        firstNode=firstNodes[blockNum]
        lastNode=lastNodes[blockNum]
        blockGraph=blockGraphs[blockNum]
        # If first mutation is covered by this primers pair, we save edge
        if amplBlockStart1<=firstNode<=amplBlockEnd1:
            blockGraph.add_edge(firstNode,primerPairName1)
            firstMut=True
            firstMutNum1=firstNode
        if amplBlockStart1<=lastNode<=amplBlockEnd1:
            blockGraph.add_edge(lastNode,primerPairName1)
            lastMut=True
            lastMutNum1=lastNode
        if lastMut and firstMut:
            blockGraph.add_edge('end',primerPairName1)                                     
            continue
        # Get index of the last mutation that is covered by this primers pair
        ## Extract the last position of amplBlock, insert it to list and get index
        if not lastMut:
            coords=chromRegionsCoords.copy()
            coords.append(amplBlockEnd1)
            lastMutNum1=sorted(coords).index(amplBlockEnd1)-1
        if not firstMut:
            coords=chromRegionsCoords.copy()
            coords.append(amplBlockStart1)
            firstMutNum1=sorted(coords).index(amplBlockStart1)
        for primerPairName2,primers2 in sorted(chromPrimersInfoByChrom.items(),key=lambda item:item[1][0][0])[i+1:]:
            if primerPairName1==primerPairName2: continue
            amplBlockStart2=primers2[0][0]+primers2[0][1]
            amplBlockEnd2=primers2[1][0]-primers2[1][1]
            # Calculate weight of this edge - distance between amplicons
            if amplBlockStart1<=amplBlockStart2:
                if lastMut: continue
                nextMut=chromRegionsCoords[lastMutNum1+1]
                if not (amplBlockStart2<=nextMut<=amplBlockEnd2):
                    continue
                # If two amplicons overlap more than user defined
                if amplBlockEnd1>=amplBlockStart2+args.maxoverlap:
                    continue
                # weight is a length of amplicons intersection or minus distance between them
                ## Less weight is a worth, the more weight is better
                weight=maxAmplLen-min(50,primers2[0][0]-primers1[1][0]) # if distance is too large (>50), leave it as 50
            else:
                if firstMut:
                    continue
                nextMut=chromRegionsCoords[max(firstMutNum1-1,0)]
                if not (amplBlockStart2<=nextMut<=amplBlockEnd2):
                    continue
                if amplBlockEnd2>=amplBlockStart1+args.maxoverlap:
                    continue
                weight=maxAmplLen-min(50,primers1[0][0]-primers2[1][0]) # if distance is too large (>50), leave it as 50
            blockGraph.add_edge(primerPairName1,
                                primerPairName2,
                                attr_dict={'weight':weight})
    blocksFinalShortestPaths=[]
    for i,blockGraph in enumerate(blockGraphs):
        if firstNodes[i]==lastNodes[i]:
            try:
                paths=tuple(nx.algorithms.shortest_paths.generic.all_shortest_paths(blockGraph,firstNodes[i],'end'))
            except KeyError:
                print('ERROR!',firstNodes[i])
                print(blockGraph.edges())
                print(chromPrimersInfoByChrom)
                print(chromRegionsCoords)
                exit(1)
            finalShortestPaths=[]
            for path in paths:
                finalShortestPaths.append(path[1:-1])
        else:
            try:
                path=tuple(nx.algorithms.shortest_paths.generic.shortest_path(blockGraph,firstNodes[i],lastNodes[i],'weight'))
            except nx.exception.NetworkXNoPath as e:
                print('ERROR! Too low value of maximal overlap (-maxoverlap) or of initially designed primers (-primernum): '+str(args.maxoverlap)+' and '+str(args.primernum1)+'. Try to increase one of them')
                logger.error(' Too low value of maximal overlap (-maxoverlap) or of initially designed primers (-primernum): '+str(args.maxoverlap)+' and '+str(args.primernum1)+'. Try to increase one of them')
                print(str(lastNodes[i])+' is not reachable from '+str(firstNodes[i]))
                logger.error(str(lastNodes[i])+' is not reachable from '+str(firstNodes[i]))
                logger.error(str(blockGraph.edges()))
                logger.error(str(chromPrimersInfoByChrom))
                exit(1)
            g1=blockGraph.subgraph(path)
            shortestPaths={path:g1.size(weight='weight')}
            minPathLen=len(path)
            maxAnalysisVars=100
            analysisNum=0
            while(True):
                analysisNum+=1
                edgelist=blockGraph.edges()
                random.shuffle(edgelist)
                g=nx.Graph(edgelist)
                path=tuple(nx.algorithms.shortest_paths.generic.shortest_path(g,firstNodes[i],lastNodes[i],'weight'))
                g1=blockGraph.subgraph(path)
                if len(path)<minPathLen:
                    shortestPaths={path:g1.size(weight='weight')}
                    minPathLen=len(path)
                    continue
                elif len(path)>minPathLen:
                    continue
                shortestPaths[path]=g1.size(weight='weight')
                if len(shortestPaths.keys())>=args.returnVariantsNum*2: break
                if analysisNum>maxAnalysisVars: break
            finalShortestPaths=[]
            j=0
            for path,value in sorted(shortestPaths.items(),key=itemgetter(1),reverse=True):
                finalShortestPaths.append(path[1:-1])
                j+=1
                if j>=args.returnVariantsNum: break
        blocksFinalShortestPaths.append(finalShortestPaths)
    return(chrom,blocksFinalShortestPaths)

def makeFinalMultiplexes(initialGraph,multiplexes=[],multNum=2):
    # We try to make cliques from it
    if len(initialGraph)==0:
        return(multiplexes)
    elif len(initialGraph.nodes())==1:
        return(multiplexes+[initialGraph.nodes()])
    elif (len(initialGraph.nodes())==multNum and
          len(multiplexes)==0):
        return([[x] for x in initialGraph.nodes()])
    graph=initialGraph.copy()
    nodesNums={}
    for node in graph.nodes():
        nodesNums[node]=graph.degree(node)
    if multNum==0:
        return(multiplexes)
    minNodesNum=math.floor(len(initialGraph)/multNum)
    maxNodesNum=math.ceil(len(initialGraph)/multNum)
    cls=[]
    edgelist=graph.edges()
    random.shuffle(edgelist)
    graph=nx.Graph(edgelist)
    for i,cl in enumerate(clique.find_cliques(graph)):
        if len(cl)>1:
            cls.append(cl)
            break
    if len(cls)==0:
        print("ERROR! Cannot choose multiplex (clique) from designed graph!")
        print('Length of input graph is:',len(initialGraph))
        print(initialGraph.edges())
        exit(1)
    if len(cls[0])>maxNodesNum:
        cls[0]=sorted(cls[0],key=graph.degree)
        cls[0]=cls[0][:-(len(cls[0])-maxNodesNum)]
    multiplexes.append(cls[0])
    graph.remove_nodes_from(cls[0])
    if len(graph)>0 and multNum-1>0:
        multiplexes=makeFinalMultiplexes(graph,multiplexes,multNum-1)
    return(multiplexes)

def sortAmpliconsToMultiplexes(globalMultiplexesContainer,globalMultiplexNums,args):
    multiplexes=[]
    leftUnsortedAmpls=[]
    for k,(key,containerNodes) in enumerate(globalMultiplexesContainer.items()):
        print('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
        logger.info('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
        mults=key.split('_')
        localMultiplexNums=nx.Graph()
        localMultiplexNums=globalMultiplexNums.subgraph(containerNodes+leftUnsortedAmpls)
        cls=[]
        cls=makeFinalMultiplexes(localMultiplexNums,[],len(mults))
        multiplexes.extend(cls)
        sumLen=sum([len(x) for x in cls])
        allSortedAmpls=[]
        for x in cls:
            allSortedAmpls.extend(x)
        if sumLen==len(localMultiplexNums):
            print(' # Number of multiplexes:',len(cls))
            logger.info(' # Number of multiplexes: '+str(len(cls)))
            for z,cl in enumerate(cls):
                print('  # Number of amplicons in multiplex',mults[z],len(cl))
                logger.info('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
        elif sum([len(x) for x in cls])<len(localMultiplexNums):
            print('Number of designed multiplexes:',len(cls))
            logger.warn(' # Number of designed multiplexes: '+str(len(cls)))
            for z,cl in enumerate(cls):
                print('Number of amplicons in multiplex',mults[z],len(cl))
                logger.warn('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
            leftNodesGraph=localMultiplexNums.copy()
            leftNodesGraph.remove_nodes_from(allSortedAmpls)
            print('  But the following '+str(len(leftNodesGraph.nodes()))+' amplicons could not be sorted to any of designed multiplex:')
            logger.warn('  But the following '+str(len(leftNodesGraph.nodes()))+' amplicons could not be sorted to any of designed multiplex:')
            for leftAmpl in leftNodesGraph.nodes():
                print(leftAmpl)
                logger.warn(leftAmpl)
            if k==len(globalMultiplexesContainer)-1:
                print('Try to change multiplex numbers in the input file.')
                logger.warn('Try to change multiplex numbers in the input file.')
            else:
                leftUnsortedAmpls=leftNodesGraph.nodes()
                print('NGS_primerplex will try to add them to the next group of multiplexes.')
                logger.warn('NGS_primerplex will try to add them to the next group of multiplexes.')
        else:
            print('UNKNOWN ERROR! Number of nodes in the final graph is more than in the initial graph!')
            print(localMultiplexNums.nodes())
            print(cls)
            logger.error('UNKNOWN ERROR! Number of nodes in the final graph is more than in the initial graph!')
            logger.error(str(localMultiplexNums.nodes()))
            logger.error(str(cls))
            exit(1)
    return(multiplexes)

# Function that checks that new primer fits one of the multiplexes
## primers is a list of two primers and their amplicon coordinates [primer1, primer2, amplStart, amplEnd]
## nums are numbers of multiplexes to which we try to put new primer
## globalMultiplexNums are lists of multiplexes with primers and their amplicon coordinates [primer1, primer2, amplStart, amplEnd]
## unspecificPrimers is a list of primer pairs that form unspecific products 
def checkPrimersFit(primers,primersToCompare,minmultdimerdg=-6,unspecificPrimers=None):
    leftPrimer,rightPrimer=primersToCompare[0:2]
    chrom,amplStart,amplEnd=primersToCompare[3:6]
    # We need to check the following:
    ## Current primers pair does not overlap with primersToCompare (the most important!)
    ## Current primers pair does not form any secondary structure (important, if with 5'-overhang)
    ## Current primers pair does not form any unspecific product
    ## Maybe, GC-content difference
    # Overlapping is the most important point of checking
    if int(chrom)==int(primers[2]):
        inter=set(range(primers[3],primers[4])).intersection(list(range(amplStart,amplEnd)))
        if len(inter)>5:
            # If there is intersection of length more than 5 bp, these primers pair do not correspond each other
            return(False,'Amplicons overlap')
    # Secondary structure with 5'-overhang
    maxdG=minmultdimerdg
    maxdG2=-12
    dG1=primer3.calcHeterodimer(primers[0],leftPrimer).dg/1000
    dG2=primer3.calcHeterodimer(primers[0],rightPrimer).dg/1000
    dG3=primer3.calcHeterodimer(primers[1],leftPrimer).dg/1000
    dG4=primer3.calcHeterodimer(primers[1],rightPrimer).dg/1000
    if ((dG1<maxdG and (revCompl(primers[0][-5:]) in leftPrimer[1:] or revCompl(leftPrimer[-5:]) in primers[0][1:]))
        or dG1<maxdG2):
        return(False,'Heterodimer of F-primer with F-primer')
    if ((dG2<maxdG and (revCompl(primers[0][-5:]) in rightPrimer[1:] or revCompl(rightPrimer[-5:]) in primers[0][1:]))
        or dG2<maxdG2):
        return(False,'Heterodimer of F-primer with R-primer')
    if ((dG3<maxdG and (revCompl(primers[1][-5:]) in leftPrimer[1:] or revCompl(leftPrimer[-5:]) in primers[1][1:]))
        or dG3<maxdG2):
        return(False,'Heterodimer of R-primer with F-primer')
    if ((dG4<maxdG and (revCompl(primers[1][-5:]) in rightPrimer[1:] or revCompl(rightPrimer[-5:]) in primers[1][1:]))
        or dG4<maxdG2):
        return(False,'Heterodimer of R-primer with R-primer')
    # Unspecific products
    if unspecificPrimers:
        if ('_'.join([primers[0],leftPrimer]) in unspecificPrimers
            or '_'.join([primers[0],rightPrimer]) in unspecificPrimers
            or '_'.join([primers[1],leftPrimer]) in unspecificPrimers
            or '_'.join([primers[1],rightPrimer]) in unspecificPrimers):
            return(False,'Unspecific product')
    return(True,None)

def extractGenomeSeq(refFa,chrom,start,end):
    attempts=0
    seq=None
    if 'chr' not in str(chrom):
        chrom='chr'+str(chrom)
    while(seq is None):
        try:
            seq=refFa.fetch(region=chrom+':'+str(start)+'-'+str(end))
        except:
            seq=None
            attempts+=1
            if attempts>=10:                    
                print('ERROR: could not extract genome sequence!',)
                logger.error('Could not extract genome sequence!')
                print(refFa.filename)
                logger.error(refFa.filename)
                print(chrom,start,end)
                logger.error('chrom:'+str(start)+'-'+str(end))
                exit(1)
    return(seq.upper())

def showPercWork(done,allWork):
    percDoneWork=round((done/allWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def revCompl(nuc):
    return(str(Seq.Seq(nuc).reverse_complement()))

# Section of input arguments
par=argparse.ArgumentParser(description='This script constructs primers for multiplex NGS panels')
par.add_argument('--regions-file','-regions',dest='regionsFile',type=str,help='file with regions for amplification in the following format:'
                 'Chromosome{Tab}Start_Position{Tab}End_Position{Tab}Amplicon_Name{Tab}\n'
                 'Desired_Multiplex_Numbers(optional){Tab}Type_Of_Primers(only left/only right/both)(optional){Tab}'
                 'Use_Whole_Region(optional)',required=True)
par.add_argument('--primers-file','-primers',dest='primersFile',type=str,help='file with previously designed internal primers. Use this parameter, if you want only to design external primers',required=False)
par.add_argument('--draft-primers','-draft',dest='draftFile',type=str,help='file with internal primers previously designed for part of input regions. The program will design primers for the left regions',required=False)
par.add_argument('--whole-genome-ref','-wgref',dest='wholeGenomeRef',type=str,help='file with INDEXED whole-genome reference sequence',required=True)
par.add_argument('--genome-version','-gv',dest='genomeVersion',type=str,help='version of genome (hg19 or hg38). Default: hg19',required=False,default='hg19')
par.add_argument('--min-amplicon-length','-minampllen',dest='minAmplLen',type=int,help='minimal length of amplicons. Default: 75',required=False,default=75)
par.add_argument('--max-amplicon-length','-maxampllen',dest='maxAmplLen',type=int,help='maximal length of amplicons. Default: 100',required=False,default=100)
par.add_argument('--optimal-amplicon-length','-optampllen',dest='optAmplLen',type=int,help='optimal length of amplicons. Default: 90',required=False,default=90)
par.add_argument('--min-primer-length','-minprimerlen',dest='minPrimerLen',type=int,help='minimal length of primers. Default: 18',required=False,default=18)
par.add_argument('--max-primer-length','-maxprimerlen',dest='maxPrimerLen',type=int,help='maximal length of primers. Default: 25',required=False,default=25)
par.add_argument('--optimal-primer-length','-optprimerlen',dest='optPrimerLen',type=int,help='optimal length of primers. Default: 23',required=False,default=23)
par.add_argument('--min-primer-melting-temp','-minprimermelt',dest='minPrimerMelt',type=int,help='minimal melting temperature of primers, degrees Celsius. Default: 62',required=False,default=62)
par.add_argument('--max-primer-melting-temp','-maxprimermelt',dest='maxPrimerMelt',type=int,help='maximal melting temperature of primers, degrees Celsius. Default: 66',required=False,default=66)
par.add_argument('--optimal-primer-melting-temp','-optprimermelt',dest='optPrimerMelt',type=int,help='optimal melting temperature of primers, degrees Celsius. Default: 64',required=False,default=64)
par.add_argument('--min-primer-gc','-minprimergc',dest='minPrimerGC',type=int,help='minimal acceptable GC-content for primers. Default: 25',required=False,default=25)
par.add_argument('--max-primer-gc','-maxprimergc',dest='maxPrimerGC',type=int,help='maximal acceptable GC-content for primers. Default: 60',required=False,default=60)
par.add_argument('--optimal-primer-gc','-optprimergc',dest='optPrimerGC',type=int,help='optimal acceptable GC-content for primers. Default: 40',required=False,default=40)
par.add_argument('--min-primer-end-gc','-minprimerendgc',dest='minPrimerEndGC',type=int,help="minimal acceptable number of G or C nucleotides within last 5 nucleotides of 3'-end of primers. Default: 1",required=False,default=1)
par.add_argument('--max-primer-end-gc','-maxprimerendgc',dest='maxPrimerEndGC',type=int,help="maximal acceptable number of G or C nucleotides within last 5 nucleotides of 3'-end of primers. Default: 3",required=False,default=3)
par.add_argument('--opt-primer-end-gc','-optprimerendgc',dest='optPrimerEndGC',type=int,help="optimal number of G or C nucleotides within last 5 nucleotides of 3'-end of primers. Default: 2",required=False,default=2)
par.add_argument('--max-primer-poly-n','-maxprimerpolyn',dest='maxPrimerPolyN',type=int,help="maximal acceptable length of some poly-N in primers. Default: 3",required=False,default=3)
par.add_argument('--max-primer-compl-end-th','-maxprimercomplendth',dest='maxPrimerComplEndTh',type=int,help="maximal Tm for complementarity of 3'-ends of primers. Default: 15",required=False,default=15)
par.add_argument('--max-primer-compl-any-th','-maxprimercomplanyth',dest='maxPrimerComplAnyTh',type=int,help="maximal Tm for any complementarity of primers. Default: 30",required=False,default=30)
par.add_argument('--max-primer-hairpin-th','-maxprimerhairpinth',dest='maxPrimerHairpinTh',type=int,help="maximal melting temperature of primer hairpin structure. Default: 40",required=False,default=40)
par.add_argument('--max-primer-nonspecific','-maxprimernonspec',dest='maxPrimerNonspec',type=int,help="maximal number of nonspecific regions to which primer can hybridizes. Default: 1000",required=False,default=1000)
par.add_argument('--max-amplicons-overlap','-maxoverlap',dest='maxoverlap',type=int,help='maximal length of overlap between two amplified blocks (it does not include primers). Default: 5',required=False,default=5)
par.add_argument('--primers-number1','-primernum1',dest='primernum1',type=int,help='number of primer that user wants to get on the 1st stage. The more this value, the more precise the choice of primers, but the longer the design time. Default: 5',required=False,default=5)
par.add_argument('--auto-adjust-parameters','-autoadjust',dest='autoAdjust',action='store_true',help='use this parameter if you want NGS-PrimerPlex to automatically use less stringent parameters if no primer were constructed for some region')
par.add_argument('--return-variants-number','-returnvariantsnum',dest='returnVariantsNum',type=int,help='number of multiplexes variants that user wants to get after all analyses and filters. Default: 1',required=False,default=1)
par.add_argument('--embedded-amplification','-embedded',dest='embeddedAmpl',action='store_true',help='use this parameter if you want to create NGS-panel with embedded amplification')
par.add_argument('--min-internal-primer-shift','-minprimershift',dest='minPrimerShift',type=int,help="minimal shift of external primer from the 3'-end of internal primer. Default: 5",required=False,default=5)
par.add_argument('--opt-external-amplicon-length','-optextampllen',dest='optExtAmplLen',type=int,help="optimal length of the external amplicons. Default: 110",required=False,default=110)
par.add_argument('--max-external-amplicon-length','-maxextampllen',dest='maxExtAmplLen',type=int,help="maximal length of the external amplicons. Default: 130",required=False,default=130)
par.add_argument('--do-blast','-blast',dest='doBlast',action='store_true',help='use this parameter if you want to perform Blast-analysis of constructed primers')
par.add_argument('--substititutions-num','-subst',dest='substNum',type=int,help='accepted number of substitutions for searching primers in genome. Default: 2',required=False,default=2)
par.add_argument('--max-nonspecific-amplicon-length','-maxnonspeclen',dest='maxNonSpecLen',type=int,help='maximal length of nonspecific amplicons that the program should consider. For example, if you design primers for DNA from serum, you can set it as 150. Default: 200',required=False,default=200)
par.add_argument('--snps','-snps',dest='snps',action='store_true',help="use this parameter if you want to check that 3'-ends of your primers do not cover any SNPs with high frequency")
par.add_argument('--min-multiplex-dimer-dg','-minmultdimerdg',dest='minMultDimerdG',type=int,help="minimal value of free energy of primer dimer formation in one multiplex in kcal/mol. Default: -6",required=False,default=-6)
par.add_argument('--snp-freq','-freq',dest='snpFreq',type=float,help='minimal frequency of SNP in whole population to consider it high-frequent SNP. Default: 0.05',required=False,default=0.05)
par.add_argument('--threads','-th',dest='threads',type=int,help='number of threads. Default: 2',required=False,default=2)
par.add_argument('--run-name','-run',dest='runName',type=str,help='name of program run. It will be used in the output file names',required=False)
par.add_argument('--skip-uncovered','-skip',dest='skipUndesigned',action='store_true',help='use this parameter if you want to skip some targets for which primers can not be designed with defined parameters')
# Parameters for calculating thermodynamic parameters
par.add_argument('--monovalent-concentration','-mv',dest='mvConc',type=int,help='Concentration of monovalent cations, commonly K+ or NH4+, in mM. Default: 50',required=False,default=50)
par.add_argument('--divalent-concentration','-dv',dest='dvConc',type=int,help='Concentration of divalent cations, commonly Mg2+, in mM. Default: 3',required=False,default=3)
par.add_argument('--dntp-concentration','-dntp',dest='dntpConc',type=float,help='Total concentration of dNTPs. If you have each dNTP with concantration 0.2 mM, then total is 0.8 mM. Default: 0.8',required=False,default=0.8)
par.add_argument('--primer-concentration','-primerconc',dest='primerConc',type=int,help='Concentration of each primer, in nM. Default: 250',required=False,default=250)
args=par.parse_args()

if not args.runName:
    args.runName=''
else:
    args.runName='_'+args.runName
runName=args.runName
inputFileBase=args.regionsFile[:-4]
# Set logging
logger=logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler=logging.FileHandler(inputFileBase+'_NGS_primerplex'+runName+'.log')
handler.setLevel(logging.INFO)
formatter=logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.info('The command was:\n'+' '.join(sys.argv))
logger.info('Arguments used:\n')
for arg in vars(args):
    logger.info(str(arg)+' '+str(getattr(args,arg)))
# Section of input arguments control
try:
    regionsFile=open(args.regionsFile)
except FileNotFoundError:
    print('ERROR! Input file was not found: '+args.regionsFile)
    logger.error('Input file was not found: '+args.regionsFile)
    exit(1)
if args.primersFile and args.draftFile:
    print('ERROR! You can use only primers file OR draft primers file. Leave one of this arguments')
    logger.error('ERROR! You can use only primers file OR draft primers file. Leave one of this arguments')
    exit(1)
inputDir=args.regionsFile[:args.regionsFile.rfind('/')+1]
wgref=args.wholeGenomeRef
if not args.maxPrimerHairpinTh:
    args.maxPrimerHairpinTh=args.minPrimerMelt-10
if (args.embeddedAmpl or args.primersFile) and args.maxExtAmplLen<=args.maxAmplLen+2*args.minPrimerShift:
    print('ERROR! Maximal length of an extrenal amplicon should be more than maximal length of internal one plus two lengths of mininimal primer shift')
    logger.error('Maximal length of an extrenal amplicon should be more than maximal length of internal one plus two lengths of mininimal primer shift')
    exit(1)
if not os.path.exists(wgref):
    print('#'*20+'\nERROR! Whole-genome reference file does not exist:',wgref)
    logger.error('Whole-genome reference file does not exist:'+wgref)
    exit(1)
if args.genomeVersion not in ['hg19','hg38']:
    print('#'*20+'\nERROR! Unaccaptable genome version was chosen:',args.genomeVersion+'. It should be hg19 or hg38')
    logger.error('WUnaccaptable genome version was chosen:'+args.genomeVersion+'. It should be hg19 or hg38')
    exit(1)

# We make primer3 parameters for each input position
## Later we will send them into multithreading pool
print('Reading input file...')
logger.info('Reading input file...')
allRegions,regionsNames,regionsCoords,regionNameToChrom,regionNameToMultiplex,regionNameToPrimerType=readInputFile(regionsFile)

# If user also used file with primers
if args.primersFile:
    print('Reading file with primers...')
    logger.info('Reading file with primers...')
    inputInternalPrimers,amplfiedRegions=readPrimersFile(args.primersFile)
    # Check that input primers cover all input regions
    chromosomesWithUncoveredRegions=[]
    for chrom,coveredRegions in amplfiedRegions.items():
        inter=coveredRegions.intersection(regionsCoords[chrom])
        if len(inter)!=len(regionsCoords[chrom]):
            chromosomesWithUncoveredRegions.append(chrom)
    if len(chromosomesWithUncoveredRegions)>0:
        print('ERROR! Some regions of the following chromosomes are not covered with input primers:')
        logger.error('Some regions of the following chromosomes are not covered with input primers:')
        for chrom in chromosomesWithUncoveredRegions:
            print(chrom)
            logger.error(chrom)
        exit(1)
    # If user wants to automatically sort primers pairs by multiplexes
    ## We create file for storing all problematic pairs of primers pairs
    if len(regionNameToMultiplex)>0:
        multiplexProblemsWB=xls.Workbook(args.primersFile[:-4]+'_amplicons_multiplex_incompatibility.xls')
        mpws=multiplexProblemsWB.add_worksheet('Internal_Primers')
        mpws.write_row(0,0,['Primers_Pair1','Primers_Pair2','Problem while joining to one multiplex'])
        colsWidth=[40,40,20]
        for k,colsWidth in enumerate(colsWidth):
            mpws.set_column(k,k,colsWidth)
    mpwsRowNum=1
    wbw=xls.Workbook(args.primersFile[:-4]+'_with_external_primers.xls')
    wsw1=wbw.add_worksheet('NGS_Primerplex_External_Primers')
    wsw1.write_row(0,0,['#','Left_Primer_Seq','Right_Primer_Seq','Amplicon_Name','Chrom','Amplicon_Start','Amplicon_End','Amplicon_Length',
                       'Amplified_Block_Start','Amplified_Block_End','Left_Primer_Tm','Right_Primer_Tm','Left_Primer_Length','Right_Primer_Length','Left_GC','Right_GC','Desired_Multiplex','Designed_Multiplex'])
    colsWidth1=[5,30,30,15,6,12,12,13,12,12,
               10,10,12,12,7,7,12,12]
    colsWidth2=[5,30,30,15,6,12,12,13,12,12,
               10,10,12,12,7,7,15,15,12,12]
    for k,colsWidth in enumerate(colsWidth1):
        wsw1.set_column(k,k,colsWidth)
    outputInternalPrimers={}
    amplNameToRowNum={}
    globalMultiplexNums=nx.Graph()
    globalMultiplexesContainer={}
    for i,internalPrimers in enumerate(inputInternalPrimers):
        amplStart=int(internalPrimers[4])
        amplEnd=int(internalPrimers[5])
        amplLen=amplEnd-amplStart+1
        amplBlockStart=int(internalPrimers[7])
        amplBlockEnd=int(internalPrimers[8])
        # We will re-compute primers Tm, because input primers parameters may have been evaluated in another conditions (salt and primers concentrations etc.)
        leftPrimerTm=int(primer3.calcTm(internalPrimers[0],mv_conc=args.mvConc,dv_conc=args.dvConc,dntp_conc=args.dntpConc,dna_conc=args.primerConc))
        rightPrimerTm=int(primer3.calcTm(internalPrimers[1],mv_conc=args.mvConc,dv_conc=args.dvConc,dntp_conc=args.dntpConc,dna_conc=args.primerConc))
        leftGC=int(internalPrimers[13])
        rightGC=int(internalPrimers[14])
        chrom=str(int(internalPrimers[3]))
        amplName=internalPrimers[2]
        regionName=amplName[:amplName.rfind('_')]
        wsw1.write_row(i+1,0,[i+1,internalPrimers[0],internalPrimers[1],amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,
                                        leftPrimerTm,rightPrimerTm,len(internalPrimers[0]),len(internalPrimers[1]),leftGC,rightGC,','.join(regionNameToMultiplex[regionName])])
        outputInternalPrimers[amplName]=[internalPrimers[0],internalPrimers[1],amplName,chrom,
                                         amplStart,amplEnd,amplLen,amplBlockStart,
                                         amplBlockEnd,leftPrimerTm,rightPrimerTm,len(internalPrimers[0]),
                                         len(internalPrimers[1]),leftGC,rightGC]
        # If user set for all regions numbers of multiplexes
        ## If it is without embedded amplification, then we do it here,
        ### but if with it, then we will do it later
        if len(regionNameToMultiplex)>0:
            amplNameToRowNum[amplName]=i+1
            num='_'.join(regionNameToMultiplex[regionName])
            if num not in globalMultiplexesContainer.keys():
                globalMultiplexesContainer[num]=[amplName]
            else:
                globalMultiplexesContainer[num].append(amplName)
            # We choose multiplex all of previously added primers that fit 
            # Get input region name that current amplicon covers
            if len(globalMultiplexNums)>0:
                for node in globalMultiplexNums.nodes():
                    try:
                        fit,problem=checkPrimersFit(internalPrimers[0:2]+[chrom,amplStart,amplEnd],outputInternalPrimers[node],args.minMultDimerdG)
                    except KeyError:
                        print('ERROR:',outputInternalPrimers.keys())
                        print(globalMultiplexNums.nodes())
                        exit(1)
                    if not fit:
                        mpws.write_row(mpwsRowNum,0,[','.join(internalPrimers[0:2]),','.join(outputInternalPrimers[node][0:2]),problem])
                        mpwsRowNum+=1
                    if fit:
                        globalMultiplexNums.add_edge(node,amplName)
                    elif amplName not in globalMultiplexNums.nodes():
                        globalMultiplexNums.add_node(amplName)
            else:
                globalMultiplexNums.add_node(amplName)
    # We use primer3, too. We extract sequences of the designed above amplicons, extend them
    ## And construct primers that will surround internal amplicons
    ### So the first step is to create primer3 input files
    print('\nCreating primer3 parameters for external primers design...')
    logger.info('Creating primer3 parameters for external primers design...')
    primer3Params,amplToChrom,amplToStartCoord=createPrimer3_parameters(allRegions,args,species='human',designedInternalPrimers=outputInternalPrimers)
    
    p=ThreadPool(args.threads)
    externalPrimersNum=0
    regionsWithoutPrimers=[]
    primersForOutput={}
    wsw2=wbw.add_worksheet('NGS_Primerplex_External_Primers')
    wsw2.write_row(0,0,['#','Left_Primer_Seq','Right_Primer_Seq','Amplicon_Name','Chrom','Amplicon_Start','Amplicon_End','Amplicon_Length',
                       'Amplified_Block_Start','Amplified_Block_End','Left_Primer_Tm','Right_Primer_Tm','Left_Primer_Length','Right_Primer_Length','Left_GC','Right_GC',"Left_3'-shift","Right_3'-shift",'Multiplex','Specificity'])
    for k,colsWidth in enumerate(colsWidth2):
        wsw2.set_column(k,k,colsWidth)
    print('Constructing external primers...')
    logger.info('Constructing external primers...')
    results=[]
    for regionName,inputParams in primer3Params.items():
        for inputParam in inputParams:
            results.append(p.apply_async(runPrimer3,(regionName,inputParam,True,args.autoAdjust)))
    outputExternalPrimers={}
    extPrimersInfo={} # This variable only for blasting designed external primers
    doneWork=0
    wholeWork=len(results)
    for res in results:
        res=res.get()
        doneWork+=1
        showPercWork(doneWork,wholeWork)
        curRegionName,primerSeqs,primersCoords,primerTms,amplLens,amplScores,primer3Params=res
        if primerSeqs==None:
            regionsWithoutPrimers.append(curRegionName)
            continue
        # leftPrimer,rightPrimer,amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,leftPrimerTm,rightPrimerTm,leftPrimerLen,rightPrimerLen,leftGC,rightGC
        # Save all found variants of external primers. Later we will choose the best ones
        for k in range(int(len(primerSeqs)/2)):
            # We do not need to create multiplexes. We only want to surround alredy designed amplicons
            ## with external primers, selecting only the best ones
            try:
                chrom=amplToChrom[curRegionName]
            except:
                print('ERROR!',amplToChrom)
                print(curRegionName)
                exit(1)
            start=amplToStartCoord[curRegionName]+primersCoords[2*k+0][0]+1
            end=amplToStartCoord[curRegionName]+primersCoords[2*k+1][0]+1
            leftGC=round(100*(primerSeqs[2*k+0].count('G')+primerSeqs[2*k+0].count('C'))/len(primerSeqs[2*k+0]),2)
            rightGC=round(100*(primerSeqs[2*k+1].count('G')+primerSeqs[2*k+1].count('C'))/len(primerSeqs[2*k+1]),2)
            if chrom=='23': chromName='X'
            elif chrom=='24': chromName='Y'
            else: chromName=str(chrom)
            out=pysam.faidx(wgref,'chr'+chromName+':'+str(start-1-100)+'-'+str(end+100))
            lines=out.split('\n')
            if lines[1]=='':
                print('ERROR! Extracted sequence has no length')
                print(chromName,start-1-100,end+100)
                exit(1)
            extendedAmplSeq=''.join(lines[1:-1]).upper()
            if curRegionName not in outputExternalPrimers.keys():
                outputExternalPrimers[curRegionName]=[[primerSeqs[2*k+0],primerSeqs[2*k+1],curRegionName+'_ext',chrom,start,end,
                                                       end-start+1,start+primersCoords[2*k][1],end-primersCoords[2*k+1][1],
                                                       primerTms[2*k+0],primerTms[2*k+1],len(primerSeqs[2*k+0]),len(primerSeqs[2*k+1]),
                                                       leftGC,rightGC,outputInternalPrimers[curRegionName][7]-(start+primersCoords[0][1])+1,
                                                       end-primersCoords[1][1]-outputInternalPrimers[curRegionName][8]+1,extendedAmplSeq]]
            else:
                outputExternalPrimers[curRegionName].append([primerSeqs[0],primerSeqs[1],curRegionName+'_ext',chrom,start,end,
                                                             end-start+1,start+primersCoords[0][1],end-primersCoords[1][1],
                                                             primerTms[0],primerTms[1],len(primerSeqs[0]),len(primerSeqs[1]),
                                                             leftGC,rightGC,outputInternalPrimers[curRegionName][7]-(start+primersCoords[0][1])+1,
                                                             end-primersCoords[1][1]-outputInternalPrimers[curRegionName][8]+1,extendedAmplSeq])
            extPrimersInfo['_'.join(primerSeqs[2*k:2*k+2])]=[chrom,start,start+primersCoords[2*k][1]-1,end-primersCoords[2*k+1][1]+1,end]
            externalPrimersNum+=1
    # Statistics of the external primers design
    if len(regionsWithoutPrimers)>0:
        regionsWithoutPrimersCounter=Counter(regionsWithoutPrimers)
        if 3 in regionsWithoutPrimersCounter.values():
            print('\n # WARNING! For',list(regionsWithoutPrimersCounter.values()).count(3),'amplicons external primers were not designed! Try less stringent parameters.')
            logger.warn(' # WARNING! For '+str(list(regionsWithoutPrimersCounter.values()).count(3))+' amplicons external primers were not designed! Try less stringent parameters.')
            for regionsWithoutPrimer,value in sorted(regionsWithoutPrimersCounter.items(),key=itemgetter(1),reverse=True):
                if value<3: break
                print('   '+regionsWithoutPrimer)
                logger.info('   '+regionsWithoutPrimer)
    print('\n # Number of designed external primers: '+str(externalPrimersNum))
    logger.info(' # Number of designed external primers: '+str(externalPrimersNum))
    p.close()
    p.join()
    # Analyzing external primers specificity        
    if args.doBlast:
        print('\nAnalyzing external primers for their specificity...')
        logger.info('Analyzing external primers for their specificity...')
        specificPrimers,primersNonSpecRegionsByChrs=checkPrimersSpecificity(args.regionsFile[:-4],extPrimersInfo,args.wholeGenomeRef,
                                                                            args.runName,args.substNum,args.threads,
                                                                            args.maxNonSpecLen,args.maxPrimerNonspec,True,str(i+1))
        print(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
        logger.info(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
    # Check external primers for covering high-frequent SNPs
    if args.snps:
        print('Analyzing external primers for covering high-frequent SNPs...')
        logger.info('Analyzing external primers for covering high-frequent SNPs...')
        p=ThreadPool(args.threads)
        results=[]
        localGenomeDB={}
        for primers,info in extPrimersInfo.items():
            primer1,primer2=primers.split('_')
            chrom=info[0]
            start1=info[1]; start2=info[3]
            end1=info[2]; end2=info[4]
            strand1=1; strand2=-1
            results.append(p.apply_async(checkPrimerForCrossingSNP,(primer1,chrom,start1,end1,strand1,localGenomeDB,args.snpFreq)))
            results.append(p.apply_async(checkPrimerForCrossingSNP,(primer2,chrom,start2,end2,strand2,localGenomeDB,args.snpFreq)))
        primersCoveringSNPs=[]
        wholeWork=len(results)
        done=0
        for res in results:
            primer,result,localGenomeDB=res.get()
            if result:
                primersCoveringSNPs.append(primer)
            done+=1
            showPercWork(done,wholeWork)
        print("\n # Number of primers covering high-frequent SNPs by 3'-end: "+str(len(primersCoveringSNPs))+'. They will be removed.')
        logger.info(" # Number of primers covering high-frequent SNPs by 3'-end: "+str(len(primersCoveringSNPs))+'. They will be removed.')
    # Now we need to remove unspecific primer pairs
    unspecificExternalPairs=0
    for curRegionName,primers in outputExternalPrimers.items():
        specificPrimersPairFound=False # Variable that says, if we have found specific primer pair among all constructed for current internal amplicon
        primersPairNotCoveringSNPs=False
        chosenVar=[]
        chrom=primers[0][3]
        for k,primer in enumerate(primers):
            start=primer[4]
            if args.doBlast and not args.snps:
                if '_'.join(primer[0:2]) in specificPrimers:
                    specificPrimersPairFound=True
                    chosenVar=k
            elif args.snps and not args.doBlast:
                if primer[0] not in primersCoveringSNPs and primer[1] not in primersCoveringSNPs:
                    primersPairNotCoveringSNPs=True
                    chosenVar=k
            elif args.snps and args.doBlast:
                if ('_'.join(primer[0:2]) in specificPrimers and
                    primer[0] not in primersCoveringSNPs and
                    primer[1] not in primersCoveringSNPs):
                    chosenVar=k
                    primersPairNotCoveringSNPs=True
                    specificPrimersPairFound=True
                    break
                # If we do not find primers that correspond to both conditions
                ## We save primers pair that specific
                elif '_'.join(primer[0:2]) in specificPrimers and not specificPrimersPairFound:
                    chosenVar=k
                    specificPrimersPairFound=True
            else:
                chosenVar=0
                break
        if args.doBlast and not specificPrimersPairFound:
            unspecificExternalPairs+=1
            chosenVar=0
        outputExternalPrimers[curRegionName]=primers[chosenVar]
        if chrom not in primersForOutput.keys():
            primersForOutput[chrom]={start:{curRegionName:primers[chosenVar]}}
        elif start not in primersForOutput[chrom].keys():
            primersForOutput[chrom][start]={curRegionName:primers[chosenVar]}
        elif curRegionName not in primersForOutput[chrom][start].keys():
           primersForOutput[chrom][start][curRegionName]=primers[chosenVar]
        else:
           primersForOutput[chrom][start][curRegionName]=primers[chosenVar]
    # If user wants to automatically sort amplicons by multiplexes
    if len(regionNameToMultiplex)>0 and args.doBlast:
        # We save primers from different pairs that form unspecific amplicons (for sorting primer pairs by multiplexs later)
        print(' Searching for nonspecific amplicons that are formed by external primers from different primer pairs...')
        logger.info(' Searching for nonspecific amplicons that are formed by external primers from different primer pairs...')
        unspecificPrimers=getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,args.maxNonSpecLen)
        print(' # Number of external primer pairs that form unspecific product:',len(unspecificPrimers))
        logger.info(' # Number of external primer pairs that form unspecific product: '+str(len(unspecificPrimers)))
    elif not args.doBlast:
        unspecificPrimers=None
    if args.doBlast and unspecificExternalPairs>0:
        print('WARNING! '+str(unspecificExternalPairs)+' external primers have nonspecific regions of hybridization or may amplify nonspecific product.')
        print('You can read this statistics in the corresponding output files.')
        logger.warn('WARNING! '+str(unspecificExternalPairs)+' external primers have nonspecific regions of hybridization or may amplify nonspecific product.')
        logger.warn('You can read this statistics in the corresponding output files.')
    # If user defined multiplex sorting
    if len(regionNameToMultiplex)>0:
        mpws=multiplexProblemsWB.add_worksheet('External_Primers')
        mpws.write_row(0,0,['Primers_Pair1','Primers_Pair2','Problem while joining to one multiplex'])
        mpwsRowNum=1
        for k,(regionName,extPrimers) in enumerate(sorted(outputExternalPrimers.items())):
            for node in sorted(outputExternalPrimers.keys())[k+1:]:
                if node==regionName: continue
                fit,problem=checkPrimersFit(extPrimers[0:2]+extPrimers[3:6],outputExternalPrimers[node],args.minMultDimerdG,unspecificPrimers)
                if not fit:
                    mpws.write_row(mpwsRowNum,0,[','.join(extPrimers[0:2]),','.join(outputExternalPrimers[node][0:2]),problem])
                    mpwsRowNum+=1
                if not fit and globalMultiplexNums.has_edge(node,regionName):
                    globalMultiplexNums.remove_edge(node,regionName)
        multiplexProblemsWB.close()
    for chrom,coords in sorted(primersForOutput.items()):
        for coord,primers in sorted(coords.items()):
            for regionName,primer in primers.items():
                try:
                    wsw2.write_row(amplNameToRowNum[regionName],0,[amplNameToRowNum[regionName]]+primer[:-1])
                    if args.doBlast and args.wholeGenomeRef:
                        if '_'.join(primer[0:2]) in specificPrimers:
                            wsw2.write(amplNameToRowNum[regionName],19,'OK')
                        else:
                            wsw2.write(amplNameToRowNum[regionName],19,'FAIL')
                except:
                    print('ERROR!',amplNameToRowNum[regionName],primer[:-1]); exit(1)
    if len(regionNameToMultiplex)>0:
        multiplexes=[]
        leftUnsortedAmpls=[]
        for k,(key,containerNodes) in enumerate(sorted(globalMultiplexesContainer.items())):
            print('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
            logger.info('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
            mults=key.split('_')
            localMultiplexNums=nx.Graph()
            localMultiplexNums=globalMultiplexNums.subgraph(containerNodes+leftUnsortedAmpls)
            cls=[]
            cls=makeFinalMultiplexes(localMultiplexNums,[],len(mults))
            multiplexes.extend(cls)
            sumLen=sum([len(x) for x in cls])
            allSortedAmpls=[]
            for x in cls:
                allSortedAmpls.extend(x)
            if sumLen==len(localMultiplexNums):
                print(' # Number of multiplexes:',len(cls))
                logger.info(' # Number of multiplexes: '+str(len(cls)))
                for z,cl in enumerate(cls):
                    print('  # Number of amplicons in multiplex',mults[z],len(cl))
                    logger.info('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
            elif sum([len(x) for x in cls])<len(localMultiplexNums):
                print('Number of designed multiplexes:',len(cls))
                logger.warn(' # Number of designed multiplexes: '+str(len(cls)))
                for z,cl in enumerate(cls):
                    print('Number of amplicons in multiplex',mults[z],len(cl))
                    logger.warn('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
                leftNodesGraph=localMultiplexNums.copy()
                leftNodesGraph.remove_nodes_from(allSortedAmpls)
                print('  But the following '+str(len(leftNodesGraph.nodes()))+' amplicons could not be sorted to any of designed multiplex:')
                logger.warn('  But the following '+str(len(leftNodesGraph.nodes()))+' amplicons could not be sorted to any of designed multiplex:')
                for leftAmpl in leftNodesGraph.nodes():
                    print(leftAmpl)
                    logger.warn(leftAmpl)
                if k==len(globalMultiplexesContainer)-1:
                    print('Try to change multiplex numbers in the input file.')
                    logger.warn('Try to change multiplex numbers in the input file.')
                else:
                    leftUnsortedAmpls=leftNodesGraph.nodes()
                    print('NGS_primerplex will try to add them to the next group of multiplexes.')
                    logger.warn('NGS_primerplex will try to add them to the next group of multiplexes.')
            else:
                print('UNKNOWN ERROR! Number of nodes in the final graph is more than in the initial graph!')
                print(localMultiplexNums.nodes())
                print(cls)
                logger.error('UNKNOWN ERROR! Number of nodes in the final graph is more than in the initial graph!')
                logger.error(str(localMultiplexNums.nodes()))
                logger.error(str(cls))
                exit(1)
        for k,multiplex in enumerate(multiplexes):
            for ampl in multiplex:
                wsw2.write(amplNameToRowNum[ampl],18,k+1)
                wsw1.write(amplNameToRowNum[ampl],17,k+1)
    print()
    logger.info('\n')
else:
    # If user use as input draft primers
    if args.draftFile:
        print('Reading file with draft primers...')
        logger.info('Reading file with draft primers...')
        # Read file with draft primers
        primersInfo,primersInfoByChrom=readDraftPrimers(args.draftFile)
        # Get regions that uncovered by already designed primers
        print('Getting positions that uncovered by draft primers...')
        logger.info('Getting positions that uncovered by draft primers...')
        uncoveredRegions,amplNames,primersToAmplNames=getRegionsUncoveredByDraftPrimers(allRegions,primersInfoByChrom)
        if len(uncoveredRegions)>0:
            # Go through all regions sorted by chromosome and coordinate of start
            print('Creating input parameters for primer3...')
            logger.info('Creating input parameters for primer3...')
            primer3Params=createPrimer3_parameters(uncoveredRegions,args,regionNameToPrimerType=regionNameToPrimerType)
            # Construct primers for each created set of parameters
            print('Constructing primers...')
            logger.info('Constructing primers...')
            # primer3Params,regionNameToChrom,args,regionsCoords=None,allRegions=None,primersInfo=None,primersInfoByChrom=None,amplNames=None,primersToAmplNames=None
            primersInfo,primersInfoByChrom,amplNames,primersToAmplNames,regionsCoords,regionNameToChrom,uncoveredRegions=constructInternalPrimers(primer3Params,regionNameToChrom,args,regionsCoords,allRegions,
                                                                                                                                                  primersInfo,primersInfoByChrom,amplNames,primersToAmplNames)
    else:                            
        # Go through all regions sorted by chromosome and coordinate of start
        print('Creating input parameters for primer3...')
        logger.info('Creating input parameters for primer3...')
        primer3Params=createPrimer3_parameters(allRegions,args,regionNameToPrimerType=regionNameToPrimerType)

        # Construct primers for each created set of parameters
        print('Constructing primers...')
        logger.info('Constructing primers...')
        primersInfo,primersInfoByChrom,amplNames,primersToAmplNames,regionsCoords,regionNameToChrom,allRegions=constructInternalPrimers(primer3Params,regionNameToChrom,args,regionsCoords,allRegions)

    # Check all primers with BWA
    if args.doBlast:
        print('Analyzing primers for their specificity...')
        logger.info('Analyzing primers for their specificity...')
        specificPrimers,primersNonSpecRegionsByChrs=checkPrimersSpecificity(args.regionsFile[:-4],primersInfo,args.wholeGenomeRef,
                                                                            args.runName,args.substNum,args.threads,
                                                                            args.maxNonSpecLen,args.maxPrimerNonspec,False,'')    
        print(' # Number of specific primer pairs:',len(specificPrimers))
        logger.info(' # Number of specific primer pairs: '+str(len(specificPrimers)))
        # Now we need to remove all unspecific primers from constructed primer pairs
        print(' Removing unspecific primer pairs...')
        logger.info(' Removing unspecific primer pairs...')
        primersInfoByChrom,amplNames=removeBadPrimerPairs(primersInfoByChrom,specificPrimers,primersToAmplNames,amplNames)
        # If we filtered primers out by specificity, we need to check that all input regions are still covered
        amplNames,allRegions,regionNameToChrom,regionsCoords=checkThatAllInputRegionsCovered(amplNames,allRegions,regionNameToChrom,regionsCoords,'by specificity',args.skipUndesigned)
        # If user wants to automatically sort amplicons by multiplexes
        if len(regionNameToMultiplex)>0:
            # We save primers from different pairs that form unspecific amplicons (for sorting primer pairs by multiplexs later)
            print(' Searching for nonspecific amplicons that are formed by primers from different primer pairs...')
            logger.info(' Searching for nonspecific amplicons that are formed by primers from different primer pairs...')
            unspecificPrimers=getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,args.maxNonSpecLen)
            print('\n # Number of primer pairs that form unspecific product:',len(unspecificPrimers))
            logger.info(' # Number of primer pairs that form unspecific product: '+str(len(unspecificPrimers)))
        writeDraftPrimers(primersInfo,args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_specificity.xls',specificPrimers)
    else:
        unspecificPrimers=None
    # Check primers for covering high-frequent SNPs
    if args.snps:
        print('Analyzing primers for covering high-frequent SNPs...')
        logger.info('Analyzing primers for covering high-frequent SNPs...')
        primerPairsNonCoveringSNPs,primersCoveringSNPs,localGenomeDB=analyzePrimersForCrossingSNP(primersInfo,args.threads,{},args.genomeVersion)
        if len(primersCoveringSNPs)>0:
            print(" Removing primer pairs covering high-frequent SNPs by 3'-end...")
            logger.info(" Removing primer pairs covering high-frequent SNPs by 3'-end...")
            primersInfoByChrom,amplNames=removeBadPrimerPairs(primersInfoByChrom,primerPairsNonCoveringSNPs,primersToAmplNames,amplNames)
            # If we filtered primers out that cover SNPs, we need to check that all input regions are still covered
            amplNames,allRegions,regionNameToChrom,regionsCoords=checkThatAllInputRegionsCovered(amplNames,allRegions,regionNameToChrom,regionsCoords,'that cover SNPs',args.skipUndesigned)
        writeDraftPrimers(primersInfo,args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_SNPs.xls',primerPairsNonCoveringSNPs)

    # Now we need to group primers into variants of multiplex PCR
    # multiplexes contains all multiplexes for ditinct joined regions
    ## [[joined_regions_1],[joined_region_2],[joined_region_3]...]
    ### Two different joined regions may be located on one or two chromsomes
    allRegionsMultiplexes={}
    pPolyN=re.compile('(A+|C+|T+|G+)')
    # We go through all chromosomes
    print('Joining primer pairs to amplified blocks...')
    logger.info('Joining primer pairs to amplified blocks...')
    p=ThreadPool(args.threads)
    results=[]
    for chrom in sorted(regionsCoords.keys()):
        results.append(p.apply_async(joinAmpliconsToBlocks,(regionsCoords[chrom],primersInfoByChrom[chrom],args.maxAmplLen,chrom)))
    doneWork=0
    wholeWork=len(results)
    showPercWork(doneWork,wholeWork)
    for res in results:
        chrom,finalShortestPaths=res.get()
        allRegionsMultiplexes[chrom]=finalShortestPaths
        doneWork+=1
        showPercWork(doneWork,wholeWork)
        
    # We need to show user statistics of amplification blocks
    totalAmplicons=0
    totalMultiplexVariants=1
    blockNum=0
    blockToChromNum={}
    # One chromosome may contain severel "amplified blocks"
    for chrom,blocks in sorted(allRegionsMultiplexes.items()):
        print('\n # Total number of amplified blocks on chromosome',str(chrom)+':',len(allRegionsMultiplexes[chrom]))
        logger.info(' # Total number of amplified blocks on chromosome '+str(chrom)+': '+str(len(allRegionsMultiplexes[chrom])))
        for i,block in enumerate(blocks):
            blockToChromNum[blockNum]=[chrom,i]
            blockNum+=1
            totalMultiplexVariants*=len(block)
            for mult in block:
                totalAmplicons+=len(mult)
    print(' # Total number of amplified blocks:',blockNum)
    print(' # Total number of amplicons:',totalAmplicons)
    print(' # Total number of multiplex variants: 10^'+str(int(round(math.log10(totalMultiplexVariants),0))))
    logger.info(' # Total number of amplified blocks: '+str(blockNum))
    logger.info(' # Total number of amplicons: '+str(totalAmplicons))
    logger.info(' # Total number of multiplex variants: 10^'+str(int(round(math.log10(totalMultiplexVariants),0))))

    # We need to select several good combinations of amplified blocks
    combinations=[]
    for i in range(args.returnVariantsNum):
        primerPairsNum=0
        # Make combination of multiplexes
        comb=[i]*blockNum
        # Go through all numbers of comb and check that there is enough multiplexes in the block
        allFit=False
        while(not allFit):
            allFit=True
            for j,c in enumerate(comb):
                chrom,chromBlockNum=blockToChromNum[j]
                if c>=len(allRegionsMultiplexes[chrom][chromBlockNum]):
                    comb[j]-=1
                    allFit=False
        combinations.append({})
        emptyMultiplexes=0
        blockNum=0
        for chrom,blocks in sorted(allRegionsMultiplexes.items()):
            combinations[-1][chrom]={}
            for block in blocks:
                multiplex=block[comb[blockNum]]
                blockNum+=1
                # Go through items of the multiplex list
                ## that is like [primers1,primers2,primers3...]
                ### each primers is primerF_primerR
                for k in range(len(multiplex)):
                    # For primersInfo[multiplex[k]] 1st [0] is a primersCoords, 2nd [0] is coordinate of F-primer; 3rd [0] is a start of F-primer
                    if primersInfo[multiplex[k]][0][0][0] not in combinations[-1][chrom].keys():
                        combinations[-1][chrom][primersInfo[multiplex[k]][0][0][0]]=multiplex[k].split('_')
                    primerPairsNum+=1
        print(' # Number of primer pairs for multiplex variant',i+1,'-',primerPairsNum)
        logger.info(' # Number of primer pairs for multiplex variant '+str(i+1)+': '+str(primerPairsNum))
    del(allRegionsMultiplexes)

    print('Writing to output...')
    logger.info('Writing to output...')
    colsWidth1=[5,30,30,15,6,12,12,13,12,12,
               10,10,12,12,7,7,12,12]
    colsWidth2=[5,30,30,15,6,12,12,13,12,12,
               10,10,12,12,7,7,15,15,12,12]
    for i,comb in enumerate(combinations):
        print('Primers combination number '+str(i+1))
        logger.info('Primers combination number '+str(i+1))
        outputInternalPrimers={}
        rFile=open(inputFileBase+'_NGS_primerplex'+runName+'_primers_combination_'+str(i+1)+'.fa','w')
        rInternalFile=open(inputFileBase+'_NGS_primerplex'+runName+'_primers_combination_'+str(i+1)+'_internal_amplicons.fa','w')
        wbw=xls.Workbook(inputFileBase+'_NGS_primerplex'+runName+'_primers_combination_'+str(i+1)+'_info.xls')
        wsw1=wbw.add_worksheet('NGS_Primerplex_Internal_Primers')
        wsw1.write_row(0,0,['#','Left_Primer_Seq','Right_Primer_Seq','Amplicon_Name','Chrom','Amplicon_Start','Amplicon_End','Amplicon_Length',
                           'Amplified_Block_Start','Amplified_Block_End','Left_Primer_Tm','Right_Primer_Tm','Left_Primer_Length','Right_Primer_Length','Left_GC','Right_GC','Desired_Multiplex','Designed_Multiplex'])
        for k,colsWidth in enumerate(colsWidth1):
            wsw1.set_column(k,k,colsWidth)
        # If user wants to automatically sort primers pairs by multiplexes
        ## We create file for storing all problematic pairs of primers pairs
        if len(regionNameToMultiplex)>0:
            multiplexProblemsWB=xls.Workbook(inputFileBase+'_NGS_primerplex'+runName+'_primers_combination_'+str(i+1)+'_amplicons_multiplex_incompatibility.xls')
            mpws=multiplexProblemsWB.add_worksheet('Internal_Primers')
            mpws.write_row(0,0,['Primers_Pair1','Primers_Pair2','Problem while joining to one multiplex'])
            colsWidth=[40,40,20]
            for k,colsWidth in enumerate(colsWidth):
                mpws.set_column(k,k,colsWidth)
            mpwsRowNum=1
        rowNum=1
        coordsNum=0
        globalMultiplexNums=nx.Graph()
        globalMultiplexesContainer={}
        amplNameToRowNum={}
        for chrom,coords in comb.items():
            for coord,c in sorted(coords.items()):
                coordsNum+=1
                rFile.write('\n'.join(['>'+primersToAmplNames['_'.join(c)][0]+'_F',
                                       c[0],
                                       '>'+primersToAmplNames['_'.join(c)][0]+'_R',
                                       c[1]])+'\n')
                if c[0]!='':
                    amplStart=primersInfo['_'.join(c)][0][0][0]
                    amplBlockStart=amplStart+primersInfo['_'.join(c)][0][0][1]
                    leftGC=round(100*(c[0].count('G')+c[0].count('C'))/len(c[0]),0)
                else:
                    amplStart=primersInfo['_'.join(c)][0][1][0]-args.maxAmplLen+1
                    amplBlockStart=amplStart+primersInfo['_'.join(c)][0][1][1]
                    leftGC=0
                if c[1]!='':
                    amplEnd=primersInfo['_'.join(c)][0][1][0]
                    amplBlockEnd=amplEnd-primersInfo['_'.join(c)][0][1][1]
                    rightGC=round(100*(c[1].count('G')+c[1].count('C'))/len(c[1]),0)
                else:
                    amplEnd=primersInfo['_'.join(c)][0][0][0]+args.maxAmplLen-1
                    amplBlockEnd=amplEnd-primersInfo['_'.join(c)][0][0][1]
                    rightGC=0
                amplLen=amplEnd-amplStart+1
                leftPrimerTm,rightPrimerTm=primersInfo['_'.join(c)][1]
                if chrom==23: chromName='X'
                elif chrom==24: chromName='Y'
                else: chromName=str(chrom) 
                out=pysam.faidx(wgref,'chr'+chromName+':'+str(amplStart-1-100)+'-'+str(amplEnd+100))
                lines=out.split('\n')
                if lines[1]=='':
                    print('ERROR! Extracted sequence has no length')
                    print(chromName,amplStart-1-100,amplEnd+100)
                    exit(1)
                extendedAmplSeq=''.join(lines[1:-1]).upper()
                rInternalFile.write('\n'.join(['>'+primersToAmplNames['_'.join(c)][0],extendedAmplSeq])+'\n')
                amplName=primersToAmplNames['_'.join(c)][0]
                regionName=amplName[:amplName.rfind('_')]
                if len(regionNameToMultiplex)>0:
                    wsw1.write_row(rowNum,0,[rowNum,c[0],c[1],amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,
                                            leftPrimerTm,rightPrimerTm,len(c[0]),len(c[1]),leftGC,rightGC,','.join(regionNameToMultiplex[regionName])])
                else:
                    wsw1.write_row(rowNum,0,[rowNum,c[0],c[1],amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,
                                            leftPrimerTm,rightPrimerTm,len(c[0]),len(c[1]),leftGC,rightGC])
                outputInternalPrimers[primersToAmplNames['_'.join(c)][0]]=[c[0],c[1],amplName,str(chrom),
                                                                           amplStart,amplEnd,amplLen,amplBlockStart,
                                                                           amplBlockEnd,leftPrimerTm,rightPrimerTm,len(c[0]),
                                                                           len(c[1]),leftGC,rightGC]
                amplNameToRowNum[amplName]=rowNum
                # If user set for all regions numbers of multiplexes
                ## If it is without embedded amplification, then we do it here,
                ### but if with it, then we will do it later
                if len(regionNameToMultiplex)>0:
                    num='_'.join(regionNameToMultiplex[regionName])
                    if num not in globalMultiplexesContainer.keys():
                        globalMultiplexesContainer[num]=[amplName]
                    else:
                        globalMultiplexesContainer[num].append(amplName)
                    # We choose multiplex all of previously added primers that fit 
                    # Get input region name that current amplicon covers
                    if len(globalMultiplexNums)>0:
                        for node in globalMultiplexNums.nodes():
                            try:
                                fit,problem=checkPrimersFit(c+[chrom,amplStart,amplEnd],outputInternalPrimers[node],args.minMultDimerdG,unspecificPrimers)
                            except KeyError:
                                print('ERROR:',outputInternalPrimers.keys())
                                print(globalMultiplexNums.nodes())
                                exit(1)
                            if not fit:
                                mpws.write_row(mpwsRowNum,0,[','.join(c),','.join(outputInternalPrimers[node][0:2]),problem])
                                mpwsRowNum+=1
                            if fit:
                                globalMultiplexNums.add_edge(node,amplName)
                            elif amplName not in globalMultiplexNums.nodes():
                                globalMultiplexNums.add_node(amplName)
                    else:
                        globalMultiplexNums.add_node(amplName)
                rowNum+=1
        rFile.close()
        rInternalFile.close()
        if len(regionNameToMultiplex)>0 and not args.embeddedAmpl:
            multiplexProblemsWB.close()
        print(' # Number of written amplicons:',coordsNum)
        logger.info(' # Number of written amplicons: '+str(coordsNum))
        # If we do not create external primers, we sort amplicons by multiplexes
        if len(regionNameToMultiplex)>0 and not args.embeddedAmpl:
            multiplexes=sortAmpliconsToMultiplexes(globalMultiplexesContainer,globalMultiplexNums,args)
            for k,multiplex in enumerate(multiplexes):
                for ampl in multiplex:
                    wsw1.write(amplNameToRowNum[ampl],17,k+1)
        
        # If we need embedded amplification
        if args.embeddedAmpl:
            # We use primer3, too. We extract sequences of the designed above amplicons, extend them
            ## And construct primers that will surround internal amplicons
            ### So the first step is to create primer3 input files
            # If user use as input draft primers
            if args.draftFile:
                print('Reading file with draft primers...')
                logger.info('Reading file with draft primers...')
                # Read file with draft primers
                extPrimersInfo,extPrimersInfoByChrom=readDraftPrimers(args.draftFile,external=True)
                # Get regions that uncovered by already designed primers
                print('Getting positions that uncovered by draft external primers...')
                logger.info('Getting positions that uncovered by draft external primers...')
                # primersInfo,primersInfoByChrom,outputInternalPrimers,args
                output=getRegionsUncoveredByDraftExternalPrimers(extPrimersInfo,extPrimersInfoByChrom,
                                                                 outputInternalPrimers,args)
                outputExternalPrimers,uncoveredInternalPrimers,amplToChrom,amplToStartCoord=output
                print('Number of internal amplicons that do not have external primers: '+str(len(uncoveredInternalPrimers)))
                logger.info('Number of internal amplicons that do not have external primers: '+str(len(uncoveredInternalPrimers)))
                if len(uncoveredInternalPrimers)>0:
                    # Go through all regions sorted by chromosome and coordinate of start
                    print('Creating input parameters for primer3...')
                    logger.info('Creating input parameters for primer3...')
                    primer3Params,amplToChrom,amplToStartCoord=createPrimer3_parameters(uncoveredRegions,args,
                                                                                        designedInternalPrimers=uncoveredInternalPrimers,
                                                                                        regionNameToPrimerType=regionNameToPrimerType)
            else:
                print('\nCreating primer3 parameters for external primers design, combination variant '+str(i+1)+'...')
                logger.info('Creating primer3 parameters for external primers design, combination variant '+str(i+1)+'...')
                primer3Params,amplToChrom,amplToStartCoord=createPrimer3_parameters(allRegions,args,species='human',
                                                                                    designedInternalPrimers=outputInternalPrimers,
                                                                                    regionNameToPrimerType=regionNameToPrimerType,
                                                                                    amplToChrom=amplToChrom,
                                                                                    amplToStartCoord=amplToStartCoord)
                outputExternalPrimers={}
                extPrimersInfo={} # This variable only for blasting designed external primers

##            createExternalPrimers(primer3Params,amplToChrom,amplToStartCoord,wgref,args,wbw,colsWidth2)
            
            p=ThreadPool(args.threads)
            externalPrimersNum=0
            regionsWithoutPrimers=[]
            primersForOutput={}
            wsw2=wbw.add_worksheet('NGS_Primerplex_External_Primers')
            wsw2.write_row(0,0,['#','Left_Primer_Seq','Right_Primer_Seq','Amplicon_Name','Chrom','Amplicon_Start','Amplicon_End','Amplicon_Length',
                               'Amplified_Block_Start','Amplified_Block_End','Left_Primer_Tm','Right_Primer_Tm','Left_Primer_Length','Right_Primer_Length','Left_GC','Right_GC',"Left_3'-shift","Right_3'-shift",'Multiplex','Specificity'])
            for k,colsWidth in enumerate(colsWidth2):
                wsw2.set_column(k,k,colsWidth)
            print('Constructing external primers...')
            logger.info('Constructing external primers...')
            results=[]
            for regionName,inputParams in primer3Params.items():
                for inputParam in inputParams:
                    results.append(p.apply_async(runPrimer3,(regionName,inputParam,True,args.autoAdjust)))
            doneWork=0
            wholeWork=len(results)
            primerDesignExplains={}
            for res in results:
                doneWork+=1
                showPercWork(doneWork,wholeWork)
                curRegionName,primerSeqs,primersCoords,primerTms,amplLens,amplScores,primer3File,designOutput=res.get()
                if primerSeqs==None:
                    regionsWithoutPrimers.append(curRegionName)
                    explanation=''
                    if 'PRIMER_PAIR_EXPLAIN' in designOutput.keys():
                        explanation=designOutput['PRIMER_PAIR_EXPLAIN']
                    elif 'PRIMER_RIGHT_EXPLAIN' in designOutput.keys():
                        explanation=designOutput['PRIMER_RIGHT_EXPLAIN']
                    elif 'PRIMER_LEFT_EXPLAIN' in designOutput.keys():
                        explanation=designOutput['PRIMER_LEFT_EXPLAIN']
                    if curRegionName not in primerDesignExplains:
                        primerDesignExplains[curRegionName]=[explanation]
                    else:
                        primerDesignExplains[curRegionName].append(explanation)
                    continue
                # leftPrimer,rightPrimer,amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,leftPrimerTm,rightPrimerTm,leftPrimerLen,rightPrimerLen,leftGC,rightGC
                # Save all found variants of external primers. Later we will choose the best ones
                for k in range(int(len(primerSeqs)/2)):
                    # We do not need to create multiplexes. We only want to surround alredy designed amplicons
                    ## with external primers, selecting only the best ones
                    try:
                        chrom=amplToChrom[curRegionName]
                    except:
                        print('ERROR!',amplToChrom)
                        print(curRegionName)
                        exit(1)
                    if primerSeqs[2*k+0]!='':
                        leftGC=round(100*(primerSeqs[2*k+0].count('G')+primerSeqs[2*k+0].count('C'))/len(primerSeqs[2*k+0]),2)
                        start=amplToStartCoord[curRegionName]+primersCoords[2*k+0][0]
                    else:
                        leftGC=0
                        start=amplToStartCoord[curRegionName]
                    if primerSeqs[2*k+1]!='':
                        rightGC=round(100*(primerSeqs[2*k+1].count('G')+primerSeqs[2*k+1].count('C'))/len(primerSeqs[2*k+1]),2)
                        end=amplToStartCoord[curRegionName]+primersCoords[2*k+1][0]
                    else:
                        rightGC=0
                        end=amplToStartCoord[curRegionName]+args.maxExtAmplLen-1
                    if chrom==23: chromName='X'
                    elif chrom==24: chromName='Y'
                    else: chromName=str(chrom)
                    out=pysam.faidx(wgref,'chr'+chromName+':'+str(start-1-100)+'-'+str(end+100))
                    lines=out.split('\n')
                    if lines[1]=='':
                        print('ERROR! Extracted sequence has no length')
                        print(chromName,start-1-100,end+100)
                        exit(1)
                    extendedAmplSeq=''.join(lines[1:-1]).upper()
                    if curRegionName not in outputExternalPrimers.keys():
                        outputExternalPrimers[curRegionName]=[[primerSeqs[2*k+0],primerSeqs[2*k+1],curRegionName+'_ext',chrom,start,end,
                                                               end-start+1,start+primersCoords[2*k][1],end-primersCoords[2*k+1][1],
                                                               primerTms[2*k+0],primerTms[2*k+1],len(primerSeqs[2*k+0]),len(primerSeqs[2*k+1]),
                                                               leftGC,rightGC,outputInternalPrimers[curRegionName][7]-(start+primersCoords[0][1])+1,
                                                               end-primersCoords[1][1]-outputInternalPrimers[curRegionName][8]+1,extendedAmplSeq]]
                    else:
                        outputExternalPrimers[curRegionName].append([primerSeqs[0],primerSeqs[1],curRegionName+'_ext',chrom,start,end,
                                                                     end-start+1,start+primersCoords[0][1],end-primersCoords[1][1],
                                                                     primerTms[0],primerTms[1],len(primerSeqs[0]),len(primerSeqs[1]),
                                                                     leftGC,rightGC,outputInternalPrimers[curRegionName][7]-(start+primersCoords[0][1])+1,
                                                                     end-primersCoords[1][1]-outputInternalPrimers[curRegionName][8]+1,extendedAmplSeq])
                    extPrimersInfo['_'.join(primerSeqs[2*k:2*k+2])]=[[[start,primersCoords[2*k][1]],[end,primersCoords[2*k+1][1]]],primerTms,end-start+1,0,chrom]
                    externalPrimersNum+=1
            # Statistics of the external primers design
            if len(regionsWithoutPrimers)>0:
                regionsWithoutPrimersCounter=Counter(regionsWithoutPrimers)
                if 3 in regionsWithoutPrimersCounter.values():
                    print(' # WARNING! For',list(regionsWithoutPrimersCounter.values()).count(3),'amplicons external primers were not designed! Try less stringent parameters.')
                    logger.warn(' # WARNING! For '+str(list(regionsWithoutPrimersCounter.values()).count(3))+' amplicons external primers were not designed! Try less stringent parameters.')
                    for regionWithoutPrimer,value in sorted(regionsWithoutPrimersCounter.items(),key=itemgetter(1),reverse=True):
                        if value<3:
                            break
                        print('   '+regionWithoutPrimer)
                        logger.info('   '+regionWithoutPrimer)
                        for explanation in primerDesignExplains[regionWithoutPrimer]:
                            print(explanation)
                            logger.info(explanation)
            print('\n # Number of designed external primers: '+str(externalPrimersNum))
            logger.info(' # Number of designed external primers: '+str(externalPrimersNum))
            p.close()
            p.join()
            # Analyzing external primers specificity        
            if args.doBlast:
                print('\nAnalyzing external primers for their specificity...')
                logger.info('Analyzing external primers for their specificity...')
                specificPrimers,primersNonSpecRegionsByChrs=checkPrimersSpecificity(args.regionsFile[:-4],extPrimersInfo,args.wholeGenomeRef,
                                                                                    args.runName,args.substNum,args.threads,
                                                                                    args.maxNonSpecLen,args.maxPrimerNonspec,True,str(i+1))
                print(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
                logger.info(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
            # Check external primers for covering high-frequent SNPs
            if args.snps:
                print('Analyzing external primers for covering high-frequent SNPs...')
                logger.info('Analyzing external primers for covering high-frequent SNPs...')
                primerPairsNonCoveringSNPs,primersCoveringSNPs,localGenomeDB=analyzePrimersForCrossingSNP(extPrimersInfo,args.threads,localGenomeDB,args.genomeVersion)
                print("\n # Number of primers covering high-frequent SNPs by 3'-end: "+str(len(primersCoveringSNPs))+'. They will be removed.')
                logger.info(" # Number of primers covering high-frequent SNPs by 3'-end: "+str(len(primersCoveringSNPs))+'. They will be removed.')
            # Now we need to remove unspecific primer pairs
            unspecificExternalPairs=0
            for curRegionName,primers in outputExternalPrimers.items():
                specificPrimersPairFound=False # Variable that says, if we have found specific primer pair among all constructed for current internal amplicon
                primersPairNotCoveringSNPs=False
                chosenVar=[]
                chrom=primers[0][3]
                for k,primer in enumerate(primers):
                    start=primer[4]
                    if args.doBlast and not args.snps:
                        if '_'.join(primer[0:2]) in specificPrimers:
                            specificPrimersPairFound=True
                            chosenVar=k
                    elif args.snps and not args.doBlast:
                        if primer[0] not in primersCoveringSNPs and primer[1] not in primersCoveringSNPs:
                            primersPairNotCoveringSNPs=True
                            chosenVar=k
                    elif args.snps and args.doBlast:
                        if ('_'.join(primer[0:2]) in specificPrimers and
                            primer[0] not in primersCoveringSNPs and
                            primer[1] not in primersCoveringSNPs):
                            chosenVar=k
                            primersPairNotCoveringSNPs=True
                            specificPrimersPairFound=True
                            break
                        # If we do not find primers that correspond to both conditions
                        ## We save primers pair that specific
                        elif '_'.join(primer[0:2]) in specificPrimers and not specificPrimersPairFound:
                            chosenVar=k
                            specificPrimersPairFound=True
                    else:
                        chosenVar=0
                        break
                if args.doBlast and not specificPrimersPairFound:
                    unspecificExternalPairs+=1
                    chosenVar=0
                outputExternalPrimers[curRegionName]=primers[chosenVar]
                if chrom not in primersForOutput.keys():
                    primersForOutput[chrom]={start:{curRegionName:primers[chosenVar]}}
                elif start not in primersForOutput[chrom].keys():
                    primersForOutput[chrom][start]={curRegionName:primers[chosenVar]}
                elif curRegionName not in primersForOutput[chrom][start].keys():
                   primersForOutput[chrom][start][curRegionName]=primers[chosenVar]
                else:
                   primersForOutput[chrom][start][curRegionName]=primers[chosenVar]
            # If user wants to automatically sort amplicons by multiplexes
            if len(regionNameToMultiplex)>0 and args.doBlast:
                # We save primers from different pairs that form unspecific amplicons (for sorting primer pairs by multiplexs later)
                print(' Searching for nonspecific amplicons that are formed by external primers from different primer pairs...')
                logger.info(' Searching for nonspecific amplicons that are formed by external primers from different primer pairs...')
                unspecificPrimers=getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,args.maxNonSpecLen)
                print(' # Number of external primer pairs that form unspecific product:',len(unspecificPrimers))
                logger.info(' # Number of external primer pairs that form unspecific product: '+str(len(unspecificPrimers)))
            elif not args.doBlast:
                unspecificPrimers=None
            if args.doBlast and unspecificExternalPairs>0:
                print('WARNING! '+str(unspecificExternalPairs)+' external primers have nonspecific regions of hybridization or may amplify nonspecific product.')
                print('You can read this statistics in the corresponding output files.')
                logger.warn('WARNING! '+str(unspecificExternalPairs)+' external primers have nonspecific regions of hybridization or may amplify nonspecific product.')
                logger.warn('You can read this statistics in the corresponding output files.')
            # If user defined multiplex sorting
            if len(regionNameToMultiplex)>0:
                mpws=multiplexProblemsWB.add_worksheet('External_Primers')
                mpws.write_row(0,0,['Primers_Pair1','Primers_Pair2','Problem while joining to one multiplex'])
                mpwsRowNum=1
                for k,(regionName,extPrimers) in enumerate(sorted(outputExternalPrimers.items())):
                    for node in sorted(outputExternalPrimers.keys())[k+1:]:
                        if node==regionName: continue
                        fit,problem=checkPrimersFit(extPrimers[0:2]+extPrimers[3:6],outputExternalPrimers[node],args.minMultDimerdG,unspecificPrimers)
                        if not fit:
                            mpws.write_row(mpwsRowNum,0,[','.join(extPrimers[0:2]),','.join(outputExternalPrimers[node][0:2]),problem])
                            mpwsRowNum+=1
                        if not fit and globalMultiplexNums.has_edge(node,regionName):
                            globalMultiplexNums.remove_edge(node,regionName)
                multiplexProblemsWB.close()
            rowNum=1
            rExternalFile=open(inputFileBase+'_NGS_primerplex'+runName+'_primers_combination_'+str(i+1)+'_external_amplicons.fa','w')
            for chrom,coords in sorted(primersForOutput.items()):
                for coord,primers in sorted(coords.items()):
                    for regionName,primer in primers.items():
                        try:
                            wsw2.write_row(amplNameToRowNum[regionName],0,[amplNameToRowNum[regionName]]+primer[:-1])
                            if args.doBlast and args.wholeGenomeRef:
                                if '_'.join(primer[0:2]) in specificPrimers:
                                    wsw2.write(amplNameToRowNum[regionName],19,'OK')
                                else:
                                    wsw2.write(amplNameToRowNum[regionName],19,'FAIL')
                            rExternalFile.write('\n'.join(['>'+primer[2],primer[-1]])+'\n')
                            rowNum+=1
                        except KeyError:
                            print('ERROR! Region name "'+regionName+'" was not found in the dictionary:')
                            print(amplNameToRowNum)
                            print(primer[:-1])
                            exit(1)
            rExternalFile.close()
            if len(regionNameToMultiplex)>0:
                multiplexes=[]
                leftUnsortedAmpls=[]
                for k,(key,containerNodes) in enumerate(sorted(globalMultiplexesContainer.items())):
                    print('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
                    logger.info('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
                    mults=key.split('_')
                    localMultiplexNums=nx.Graph()
                    localMultiplexNums=globalMultiplexNums.subgraph(containerNodes+leftUnsortedAmpls)
                    cls=[]
                    cls=makeFinalMultiplexes(localMultiplexNums,[],len(mults))
                    multiplexes.extend(cls)
                    sumLen=sum([len(x) for x in cls])
                    allSortedAmpls=[]
                    for x in cls:
                        allSortedAmpls.extend(x)
                    if sumLen==len(localMultiplexNums):
                        print(' # Number of multiplexes:',len(cls))
                        logger.info(' # Number of multiplexes: '+str(len(cls)))
                        for z,cl in enumerate(cls):
                            print('  # Number of amplicons in multiplex',mults[z],len(cl))
                            logger.info('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
                    elif sum([len(x) for x in cls])<len(localMultiplexNums):
                        print('Number of designed multiplexes:',len(cls))
                        logger.warn(' # Number of designed multiplexes: '+str(len(cls)))
                        for z,cl in enumerate(cls):
                            print('Number of amplicons in multiplex',mults[z],len(cl))
                            logger.warn('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
                        leftNodesGraph=localMultiplexNums.copy()
                        leftNodesGraph.remove_nodes_from(allSortedAmpls)
                        print('  But the following '+str(len(leftNodesGraph.nodes()))+' amplicons could not be sorted to any of designed multiplex:')
                        logger.warn('  But the following '+str(len(leftNodesGraph.nodes()))+' amplicons could not be sorted to any of designed multiplex:')
                        for leftAmpl in leftNodesGraph.nodes():
                            print(leftAmpl)
                            logger.warn(leftAmpl)
                        if k==len(globalMultiplexesContainer)-1:
                            print('Try to change multiplex numbers in the input file.')
                            logger.warn('Try to change multiplex numbers in the input file.')
                        else:
                            leftUnsortedAmpls=leftNodesGraph.nodes()
                            print('NGS_primerplex will try to add them to the next group of multiplexes.')
                            logger.warn('NGS_primerplex will try to add them to the next group of multiplexes.')
                    else:
                        print('UNKNOWN ERROR! Number of nodes in the final graph is more than in the initial graph!')
                        print(localMultiplexNums.nodes())
                        print(cls)
                        logger.error('UNKNOWN ERROR! Number of nodes in the final graph is more than in the initial graph!')
                        logger.error(str(localMultiplexNums.nodes()))
                        logger.error(str(cls))
                        exit(1)
                for k,multiplex in enumerate(multiplexes):
                    for ampl in multiplex:
                        wsw2.write(amplNameToRowNum[ampl],18,k+1)
                        wsw1.write(amplNameToRowNum[ampl],17,k+1)
        wbw.close()
        print()
        logger.info('\n')
##wbw.close()


# TODO:
        
