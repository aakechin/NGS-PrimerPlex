#!/usr/bin/python3
# This script constructs primers for multiplex NGS panels

import argparse
import re
import os
import sys
import math
import random
import primer3
import logging
import pysam
import xlrd
import numpy
from copy import deepcopy
from multiprocessing.pool import ThreadPool
from Bio import SeqIO,Seq,pairwise2
import subprocess as sp
from operator import itemgetter
import xlsxwriter as xls
import networkx as nx
import networkx.algorithms.clique as clique
import networkx.algorithms.shortest_paths.weighted as weighted_shortest_paths
from itertools import islice
from collections import Counter

global thisDir,nameToNum,numToName
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

# Section of functions
def chrToChr(args,refFa):
    global nameToNum,numToName
    nameToNum={}
    numToName={}
##    refFa=pysam.FastaFile(args.wholeGenomeRef)
    for i,ch in enumerate(refFa.references):
        nameToNum[ch]=i+1
        numToName[i+1]=ch
    return(nameToNum,numToName)

def readInputFile(regionsFile):
    global nameToNum,numToName
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
        # chrom is a name of chromosome (e.g. chr1, chr2 etc)
        # chromInt is a number of chromosome in the reference genome file (e.g. 1,2,3 etc)
        # But the last one can be different from names (e.g. chr1 can be chromosome #2)
        # So the user shouldn't see this numbers in order not to shock him or her
        chrom=cols[0]
        try:
            chromInt=int(chrom)
        except ValueError:
            if chrom in nameToNum.keys():
                chromInt=nameToNum[chrom]
            else:
                print('ERROR (56)! Incorrect format of chromosome!')
                print('According to the defined reference genome, chromosome names are:')
                print(sorted(nameToNum.keys()))
                print('Or you can use numbers of chromosome, and they correspond to the following names:')
                for num,name in sorted(numToName.items(),
                                       key=lambda item:item[0]):
                    print(num,name)
                exit(56)
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
                if chromInt not in regionsCoords.keys() or i not in regionsCoords[chromInt]:
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
                            exit(2)
                        regionNameToPrimerType[curRegionName]=cols[5]
                    if len(cols)>6 and cols[6]!='':
                        if cols[6]=='W':
                            endShift=regEnd-regStart
                            print('ATTENTION! For region',regionName,'whole region was chosen!')
                            logger.info('ATTENTION! For region '+regionName+' whole region was chosen!')
                        else:
                            endShift=0
                    else:
                        endShift=0
                    if chrom not in allRegions.keys():
                        allRegions[chrom]={curRegionName:[chrom,i,i+endShift,curRegionName]}
                        regionsCoords[chromInt]=[i]
                        uniquePointRegions+=1
                    else:
                        # Addtionally check that all input regions are unique
                        if [chrom,i,i+endShift,curRegionName] not in allRegions[chrom].values():
                            allRegions[chrom][curRegionName]=[chrom,i,i+endShift,curRegionName]
                            regionsCoords[chromInt].append(i)
                            uniquePointRegions+=1
                    if endShift>0:
                        for i in range(regStart+1,regEnd+1):
                            regionsCoords[chromInt].append(i)
                        break
            except ValueError:
                print('ERROR: Incorrect format of input file!')
                print('It should have the following format:')
                print('Chromosome{Tab}Start_Position{Tab}End_Position{Tab}Amplicon_Name{Tab}\n'
                      'Desired_Multiplex_Numbers(optional){Tab}Type_Of_Primers(only left/only right/both)(optional){Tab}'
                      'Use_Whole_Region(optional)')
                print('But your file have the following format:')
                print(string)
                exit(4)
    # Sort dict with regions coords
    for chromInt in regionsCoords.keys():
        regionsCoords[chromInt]=sorted(regionsCoords[chromInt])
    print(' # Total number of input point regions:',totalPointInputRegions)
    print(' # Number of unique input point regions:',uniquePointRegions)
    logger.info(' # Total number of input point regions: '+str(totalPointInputRegions))
    logger.info(' # Number of unique input point regions: '+str(uniquePointRegions))
    # allRegions has chromosome names
    # regionsCoords has chromosome numbers
    # regionNameToChrom has chromosome names
    return(allRegions,regionsNames,regionsCoords,regionNameToChrom,regionNameToMultiplex,regionNameToPrimerType)

def readPrimersFile(primersFile):
    wb=xlrd.open_workbook(primersFile)
    ws=wb.sheet_by_index(0)
    internalPrimers=[]
    amplifiedRegions={}
    for i in range(ws.nrows):
        row=ws.row_values(i)
        if row[0]=='':
            break
        if i==0:
            # Check that input file has necessary format
            if row[:16]!='#	Left_Primer_Seq	Right_Primer_Seq	Amplicon_Name	Chrom	Amplicon_Start	Amplicon_End	Amplicon_Length	Amplified_Block_Start	Amplified_Block_End	Left_Primer_Tm	Right_Primer_Tm	Left_Primer_Length	Right_Primer_Length	Left_GC	Right_GC'.split('\t'):
                print('ERROR! Input file with primers has incorrect format. You can use only files with primers that has format of NGS-primerplex')
                logger.error('Input file with primers has incorrect format. You can use only files with primers that has format of NGS-primerplex')
                exit(5)
            continue
        # Make chromosome from excel cell
        chrom=getChrNum(row[4])
        if chrom not in amplifiedRegions.keys():
            amplifiedRegions[chrom]=set(range(int(row[8]),int(row[9])+1))
        else:
            amplifiedRegions[chrom].update(set(range(int(row[8]),int(row[9])+1)))
        internalPrimers.append(row[1:16])
    return(internalPrimers,amplifiedRegions)

def createPrimer3_parameters(pointRegions,args,refFa,
                             designedInternalPrimers=None,
                             regionNameToPrimerType=None,
                             regionNameNeedToBeWholeLen=None,
                             regionNameToChrom={},
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
    primerTags=deepcopy(templatePrimerTags)
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
    if designedInternalPrimers==None:
        if args.leftAdapter:
            primerTags['PRIMER_LEFT_ADAPTER']=args.leftAdapter.upper()
        if args.rightAdapter:
            primerTags['PRIMER_RIGHT_ADAPTER']=args.rightAdapter.upper()
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
                    print('ERROR #6: Unknown type of primers is necessary to be designed for the following region:')
                    logger.error('#6 Unknown type of primers is necessary to be designed for the following region:')
                    print(amplName)
                    logger.error(amplName)
                    print('This value can be only "L" (only left primer), "R" (only right primer), "B" (both primers) or nothing (both primers)')
                    logger.error('This value can be only "L" (only left primer), "R" (only right primer), "B" (both primers) or nothing (both primers)')
                    exit(6)
            regionNameToChrom[amplName]=chrom
            primer3Params[amplName]=[]
            primerTags['PRIMER_PRODUCT_SIZE_RANGE']=[[amplLen+2*args.minPrimerShift,args.maxExtAmplLen]]
            primerTags['PRIMER_PRODUCT_OPT_SIZE']=args.optExtAmplLen
            # We can shift our external amplicon on the (maximal length of extAmpl) - (amplLen)-(minPrimerShift)
            chrTargetSeqStart=amplBlockEnd-(args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen)
            # In COSMIC database chromosome X and Y are designated as 23 and 24, respectively
            ## So we need to check if primers are designed for human genome and chromosomes may be 23 and 24 instead of X and Y
            chromInt=nameToNum[chrom]
            chromName=chrom
            if (chrTargetSeqStart<0 or
                amplBlockEnd+args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen<0):
                print('ERROR #54: Unknown error with coordinates for sequence extraction from genome:')
                logger.error('#54 Unknown error with coordinates for sequence extraction from genome:')
                print(chromName,chrTargetSeqStart,amplBlockEnd,
                      args.maxExtAmplLen,args.minPrimerShift,
                      args.minPrimerLen)
                logger.error(', '.join(map(str,[chromName,chrTargetSeqStart,amplBlockEnd,
                                                     args.maxExtAmplLen,args.minPrimerShift,
                                                     args.minPrimerLen])))
                exit(54)
            regionSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                       chromName,chrTargetSeqStart,
                                       amplBlockEnd+(args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen)+1)
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
                primer3Params[amplName].append([deepcopy(seqTags),deepcopy(primerTags)])
        return(primer3Params,regionNameToChrom,amplToStartCoord)
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
                        exit(8)
                # We do not need to split regions onto several blocks
                ## because previously we splited it onto point positions
                primer3Params[curRegionName]=[]
                # In COSMIC database chromosome X and Y are designated as 23 and 24, respectively
                ## So we need to check if primers are designed for human genome and chromosomes may be 23 and 24 instead of X and Y
                chromInt=nameToNum[chrom]
                chromName=chrom
                seqTags={}
                seqTags['SEQUENCE_ID']='NGS_primerplex_'+curRegionName
                if (start-args.maxAmplLen<0 and
                    end+args.maxAmplLen>refFa.lengths[refFa.references.index(chromName)]):
                    print('ERROR (53)! The emplicon length is more then chromosome length')
                    logger.error('(53) The emplicon length is more then chromosome length')
                    print(chromName,start,args.maxAmplLen,end)
                    logger.error(', '.join([chromName,str(start),str(args.maxAmplLen),str(end)]))
                elif start-args.maxAmplLen<0:
                    regionSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                               chromName,1,end+args.maxAmplLen)
                    seqTags['SEQUENCE_TEMPLATE']=str(regionSeq)
                    prevEnd=end
                    primerPairOkRegion=[[0,start-1,
                                         end+1,args.maxAmplLen-1]]
                    seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                    primer3Params[curRegionName].append([deepcopy(seqTags),deepcopy(primerTags)])
                elif end+args.maxAmplLen>refFa.lengths[refFa.references.index(chromName)]:
                    regionSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                               chromName,max(1,start-args.maxAmplLen),
                                               refFa.lengths[refFa.references.index(chromName)])
                    seqTags['SEQUENCE_TEMPLATE']=str(regionSeq)
                    prevEnd=end
                    primerPairOkRegion=[[0,args.maxAmplLen-args.maxPrimerLen-2,
                                         args.maxAmplLen+1+end-start,
                                         len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                    seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                    primer3Params[curRegionName].append([deepcopy(seqTags),deepcopy(primerTags)])
                else:
                    regionSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                               chromName,max(1,start-args.maxAmplLen),end+args.maxAmplLen)
                    seqTags['SEQUENCE_TEMPLATE']=str(regionSeq)
                    if start==prevEnd+1:                    
                        primerPairOkRegion=[[args.maxAmplLen-args.maxPrimerLen-2,args.maxPrimerLen+2,
                                             args.maxAmplLen+1+end-start,len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                        seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                        primer3Params[curRegionName].append([deepcopy(seqTags),deepcopy(primerTags)])
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
                                primerPairOkRegion=[[args.maxAmplLen-args.maxPrimerLen-2,args.maxPrimerLen+2,
                                                     args.maxAmplLen+1+end-start,
                                                     len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                            else:
                                primerPairOkRegion=[[0,args.maxAmplLen-1,
                                                     args.maxAmplLen+1+end-start,
                                                     len(regionSeq)-args.maxAmplLen-1-(end-start)]]
                            seqTags['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']=primerPairOkRegion
                            primer3Params[curRegionName].append([deepcopy(seqTags),deepcopy(primerTags)])
                    prevEnd=end
        return(primer3Params)

def constructInternalPrimers(primer3Params,regionNameToChrom,
                             args,regionsCoords=None,
                             allRegions=None,primersInfo=None,
                             primersInfoByChrom=None,
                             amplNames=None,primersToAmplNames=None):
    # chrom is string
    p=ThreadPool(args.threads)
    # Dictionary for storing primers' info
    # Contains chromosome names
    if primersInfo==None:
        primersInfo={}
    # Dictionary for storing primers' info but primers are splitted by chromosome location
    # Contains chromosome numbers
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
            results.append(p.apply_async(runPrimer3,(regionName,inputParam,False,args)))
        showPercWork(i+1,wholeWork,args.gui)
    print()
    doneWork=0
    wholeWork=len(results)
    primerDesignExplains={}
    internalExplainToWords=['left min end GC',
                            'right min end GC',
                            'left hairpin end3',
                            'right hairpin end3',
                            'left homodimer',
                            'right homodimer',
                            'heterodimer',
                            'left homodimer end3',
                            'right homodimer end3',
                            'heterodimer end3']
    for res in results:
        doneWork+=1
        showPercWork(doneWork,wholeWork,args.gui)
        curRegionName,primerSeqs,primersCoords,primerTms,amplLens,amplScores,primer3File,designOutput=res.get()
        if curRegionName not in amplNames.keys():
            amplNames[curRegionName]=[]
        if primerSeqs==None:
            explanation=''
            if 'PRIMER_PAIR_EXPLAIN' in designOutput.keys():
                explanation=designOutput['PRIMER_PAIR_EXPLAIN']
                if explanation=='considered 0, ok 0':
                    explanation='considered 0 pairs, ok 0;'
                    if 'PRIMER_LEFT_EXPLAIN' in designOutput.keys():
                        if 'ok 0' in designOutput['PRIMER_LEFT_EXPLAIN']:
                            explanation+='\n  Left: '+designOutput['PRIMER_LEFT_EXPLAIN']
                    if 'PRIMER_RIGHT_EXPLAIN' in designOutput.keys():
                        if 'ok 0' in designOutput['PRIMER_RIGHT_EXPLAIN']:
                            explanation+='\n  Right: '+designOutput['PRIMER_RIGHT_EXPLAIN']
            elif 'PRIMER_RIGHT_EXPLAIN' in designOutput.keys():
                explanation=designOutput['PRIMER_RIGHT_EXPLAIN']
            elif 'PRIMER_LEFT_EXPLAIN' in designOutput.keys():
                explanation=designOutput['PRIMER_LEFT_EXPLAIN']
            if 'INTERNAL_EXPLAIN' in designOutput.keys():
                for k,filteredNum in enumerate(designOutput['INTERNAL_EXPLAIN']):
                    if filteredNum>0:
                        explanation+=' \n'+str(filteredNum)+' pairs were filtered by '+internalExplainToWords[k]
            if curRegionName not in primerDesignExplains:
                primerDesignExplains[curRegionName]=[]
            primerDesignExplains[curRegionName].append(explanation)
            continue
        chromName=regionNameToChrom[curRegionName]
        chrom=chromName
        chromInt=nameToNum[chromName]
        if chromInt not in primersInfoByChrom.keys():
            primersInfoByChrom[chromInt]={}
        # Extract start and end of target region
        try:
            start,end=allRegions[chromName][curRegionName][1:3]
        except KeyError:
            print('ERROR!',allRegions.keys())
            exit(10)
        # Go through each pair of primers and save them and info about them
        for i in range(int(len(primerSeqs)/2)):
            totalPrimersNum+=i+1
            # If this pair of primers has been already processed, then we save all possible info about it
            ## Including information that this primers cover other input positions
            if '_'.join(primerSeqs[2*i:2*i+2]) in primersInfo.keys():
                continue
            # Calculating coordinates of primers on a chromosome
            # If it is near the start of a chromosome
            if start-args.maxAmplLen<0:
                primersCoords[2*i][0]=primersCoords[2*i][0]
                ## If we constructed only left primer
                if primersCoords[2*i+1][0]==0:
                    # We set coordinate of the right primer as the most right position of this amplicon
                    primersCoords[2*i+1][0]=args.maxAmplLen
                else:
                    primersCoords[2*i+1][0]=primersCoords[2*i+1][0]
            else:
                primersCoords[2*i][0]=primersCoords[2*i][0]+start-args.maxAmplLen # We do not substract 1, because we want to get real coordinate, not number of symbol in coordinate
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
            if '_'.join(primerSeqs[2*i:2*i+2]) not in primersInfoByChrom[chromInt].keys():
                primersInfoByChrom[chromInt]['_'.join(primerSeqs[2*i:2*i+2])]=primersCoords[2*i:2*i+2]
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
        for regionWithoutPrimer in sorted(regionsWithoutPrimers,
                                          key=lambda key:splitNameForSorting(key)):
            if args.skipUndesigned:
                regionsCoords[nameToNum[regionNameToChrom[regionWithoutPrimer]]].remove(allRegions[regionNameToChrom[regionWithoutPrimer]][regionWithoutPrimer][1])
                if len(regionsCoords[nameToNum[regionNameToChrom[regionWithoutPrimer]]])==0:
                    regionsCoords.pop(nameToNum[regionNameToChrom[regionWithoutPrimer]])
                allRegions[regionNameToChrom[regionWithoutPrimer]].pop(regionWithoutPrimer)
                if len(allRegions[regionNameToChrom[regionWithoutPrimer]])==0:
                    allRegions.pop(regionNameToChrom[regionWithoutPrimer])
                regionNameToChrom.pop(regionWithoutPrimer)
                amplNames.pop(regionWithoutPrimer)
            print('   '+regionWithoutPrimer)
            logger.info('   '+regionWithoutPrimer)
            print(primer3Params[regionWithoutPrimer][0][0]['SEQUENCE_TEMPLATE'])
            logger.info(primer3Params[regionWithoutPrimer][0][0]['SEQUENCE_TEMPLATE'])
            if regionWithoutPrimer not in primerDesignExplains.keys():
                print('All primer pairs were filtered out by non-primer3 parameters. '
                      "Try changing NGS-PrimerPlex parameters (e.g. homodimer stability with hybridized 3'-end)")
                logger.info('All primer pairs were filtered out by non-primer3 parameters. '
                            "Try changing NGS-PrimerPlex parameters (e.g. homodimer stability with hybridized 3'-end)")
                continue
            for explanation in primerDesignExplains[regionWithoutPrimer]:
                print(explanation)
                logger.info(explanation)
        if not args.skipUndesigned:
            print(' You should use less stringent parameters')
            logger.info(' You should use less stringent parameters')
            exit(12)
    print(' # Total number of constructed primers:',totalPrimersNum)
    print(' # Total number of different constructed primers:',totalDifPrimersNum)
    logger.info(' # Total number of constructed primers: '+str(totalPrimersNum))
    logger.info(' # Total number of different constructed primers: '+str(totalDifPrimersNum))
    return(primersInfo,primersInfoByChrom,amplNames,primersToAmplNames,regionsCoords,regionNameToChrom,allRegions)    

def runPrimer3(regionName,inputParams,extPrimer,args):
    autoAdjust=args.autoAdjust
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
            exit(13)
        numReturned=out['PRIMER_PAIR_NUM_RETURNED']
##        print(numReturned)
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
            while(i<int(len(primers)/2) and len(primers)>0):
                leftEnd3_rc=str(Seq.Seq(primers[2*i][-4:]).reverse_complement())
                rightEnd3_rc=str(Seq.Seq(primers[2*i+1][-4:]).reverse_complement())
                if 'PRIMER_LEFT_ADAPTER' in primerTags.keys():
                    leftPrimer=primerTags['PRIMER_LEFT_ADAPTER']+primers[2*i]
                else:
                    leftPrimer=primers[2*i]
                if 'PRIMER_RIGHT_ADAPTER' in primerTags.keys():
                    rightPrimer=primerTags['PRIMER_RIGHT_ADAPTER']+primers[2*i+1]
                else:
                    rightPrimer=primers[2*i]
                leftHairpin=primer3.calcHairpin(leftPrimer,
                                                mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                rightHairpin=primer3.calcHairpin(rightPrimer,
                                                 mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                 dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                 dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                 dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                leftHomodimer=primer3.calcHomodimer(leftPrimer,
                                                    mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                    dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                    dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                    dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                rightHomodimer=primer3.calcHomodimer(rightPrimer,
                                                     mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                     dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                     dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                     dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                heterodimer=primer3.calcHeterodimer(leftPrimer,rightPrimer,
                                                    mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                    dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                    dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                    dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                leftHairpinEnd3=calcThreeStrikeEndHairpin(leftPrimer,
                                                          mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                          dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                          dntp_conc=primerTags['PRIMER_DNTP_CONC'],
                                                          dna_conc=primerTags['PRIMER_DNA_CONC'])
                rightHairpinEnd3=calcThreeStrikeEndHairpin(rightPrimer,
                                                          mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                          dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                          dntp_conc=primerTags['PRIMER_DNTP_CONC'],
                                                          dna_conc=primerTags['PRIMER_DNA_CONC'])
                leftHomodimerEnd3=calcThreeStrikeEndDimer(leftPrimer,leftPrimer,
                                                          mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                          dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                          dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                          dntp_conc=primerTags['PRIMER_DNTP_CONC'])
                rightHomodimerEnd3=calcThreeStrikeEndDimer(rightPrimer,rightPrimer,
                                                           mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                           dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                           dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                           dntp_conc=primerTags['PRIMER_DNTP_CONC'])
                heterodimerEnd3=calcThreeStrikeEndDimer(leftPrimer,rightPrimer,
                                                        mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                        dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                        dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                        dntp_conc=primerTags['PRIMER_DNTP_CONC'])
                internalChecks=[primers[2*i][-5:].count('G')+primers[2*i][-5:].count('C')<minEndGC,
                                primers[2*i+1][-5:].count('G')+primers[2*i+1][-5:].count('C')<minEndGC,
                                leftHairpinEnd3<-2,
                                rightHairpinEnd3<-2,
                                leftHomodimer<args.minMultDimerdG2,
                                rightHomodimer<args.minMultDimerdG2,
                                heterodimer<args.minMultDimerdG2,
                                leftHomodimerEnd3<args.minMultDimerdG1,
                                rightHomodimerEnd3<args.minMultDimerdG1,
                                heterodimerEnd3<args.minMultDimerdG1]
                if True in internalChecks:
                    out['INTERNAL_EXPLAIN']=[0]*10
                    for k,check in enumerate(internalChecks):
                        if check:
                            out['INTERNAL_EXPLAIN'][k]+=1
                            break
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
            while(i<int(len(primers)/2) and len(primers)>0):
                if 'PRIMER_LEFT_ADAPTER' in primerTags.keys():
                    leftPrimer=primerTags['PRIMER_LEFT_ADAPTER']+primers[2*i]
                else:
                    leftPrimer=primers[2*i]
                leftEnd3_rc=str(Seq.Seq(primers[2*i][-4:]).reverse_complement())
                leftHairpin=primer3.calcHairpin(leftPrimer,
                                                mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                leftHomodimer=primer3.calcHomodimer(leftPrimer,
                                                    mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                    dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                    dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                    dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                leftHairpinEnd3=calcThreeStrikeEndHairpin(leftPrimer,
                                                          mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                          dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                          dntp_conc=primerTags['PRIMER_DNTP_CONC'],
                                                          dna_conc=primerTags['PRIMER_DNA_CONC'])
                leftHomodimerEnd3=calcThreeStrikeEndDimer(leftPrimer,leftPrimer,
                                                          mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                          dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                          dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                          dntp_conc=primerTags['PRIMER_DNTP_CONC'])
                if (primers[2*i][-5:].count('G')+primers[2*i][-5:].count('C')<minEndGC
                    or leftHairpinEnd3<-2
                    or leftHomodimer<args.minMultDimerdG2
                    or leftHomodimerEnd3<args.minMultDimerdG1):
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
            while(i<int(len(primers)/2) and len(primers)>0):
                if 'PRIMER_RIGHT_ADAPTER' in primerTags.keys():
                    rightPrimer=primerTags['PRIMER_RIGHT_ADAPTER']+primers[2*i+1]
                else:
                    rightPrimer=primers[2*i]
                rightEnd3_rc=str(Seq.Seq(primers[2*i+1][-4:]).reverse_complement())
                rightHairpin=primer3.calcHairpin(rightPrimer,
                                                 mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                 dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                 dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                 dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                rightHomodimer=primer3.calcHomodimer(rightPrimer,
                                                     mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                     dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                     dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                     dntp_conc=primerTags['PRIMER_DNTP_CONC']).dg/1000
                rightHairpinEnd3=calcThreeStrikeEndHairpin(rightPrimer,
                                                          mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                          dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                          dntp_conc=primerTags['PRIMER_DNTP_CONC'],
                                                          dna_conc=primerTags['PRIMER_DNA_CONC'])
                rightHomodimerEnd3=calcThreeStrikeEndDimer(rightPrimer,rightPrimer,
                                                           mv_conc=primerTags['PRIMER_SALT_MONOVALENT'],
                                                           dv_conc=primerTags['PRIMER_SALT_DIVALENT'],
                                                           dna_conc=primerTags['PRIMER_DNA_CONC'],
                                                           dntp_conc=primerTags['PRIMER_DNTP_CONC'])
                if (primers[2*i+1][-5:].count('G')+primers[2*i+1][-5:].count('C')<minEndGC
                    or rightHairpinEnd3<-2
                    or rightHomodimer<args.minMultDimerdG2
                    or rightHomodimerEnd3<args.minMultDimerdG1):
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
##                print(primerTags)
##                print(seqTags)
##                print(sorted(seqTags.items()))
##                print(out)
##                input()
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
    try:
        wb=xlrd.open_workbook(draftFile)
    except FileNotFoundError:
        print('ERROR (46)! The defined file with draft primers was not found')
        print(draftFile)
        exit(46)
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
        chrom=getChrNum(row[9])
        chromInt=nameToNum[chrom]
        primersInfo[row[0]]=[[[int(row[1]),int(row[2])],[int(row[3]),int(row[4])]],
                             [int(row[5]),int(row[6])],
                             int(row[7]),int(row[8]),str(chrom)]
        if chromInt not in primersInfoByChrom.keys():
            primersInfoByChrom[chromInt]={row[0]:[[int(row[1]),int(row[2])],[int(row[3]),int(row[4])]]}
        elif row[0] not in primersInfoByChrom[chromInt].keys():
            primersInfoByChrom[chromInt][row[0]]=[[int(row[1]),int(row[2])],[int(row[3]),int(row[4])]]
    return(primersInfo,primersInfoByChrom)

def getChrNum(chrom):
    try:
        chrom=str(int(round(float(chrom),0)))
    except ValueError:
        if chrom not in nameToNum.keys():
            print('ERROR (55)! Chromosome '+ chrom+' not found in reference file')
            exit(55)
    return(chrom)

def getRegionsUncoveredByDraftPrimers(allRegions,primersInfoByChrom):
    # allRegions[chrom]={curRegionName:[chrom,i,i+endShift,curRegionName]}
    # primersInfoByChrom[int(chrom)]['_'.join(primerSeqs[2*i:2*i+2])]=primersCoords[2*i:2*i+2]
    amplNames={}
    primersToAmplNames={}
    uncoveredRegions=deepcopy(allRegions)
    for chromInt,coords in primersInfoByChrom.items():
        coordToRegionName={}
        chrom=numToName[chromInt]
        # If current primer pair does not cover any of target regions
        if chrom not in allRegions.keys():
            continue
        for regionCoords in allRegions[chrom].values():
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
                    if (curRegionName in uncoveredRegions[chrom].keys() and
                        allRegions[chrom][curRegionName][1]==allRegions[chrom][curRegionName][2]):
                        uncoveredRegions[chrom].pop(curRegionName)
                    elif (curRegionName in uncoveredRegions[chrom].keys() and
                          coord[0][0]+coord[0][1]<=allRegions[chrom][curRegionName][1] and
                          coord[1][0]-coord[1][1]>=allRegions[chrom][curRegionName][2]):
                        uncoveredRegions[chrom].pop(curRegionName)
                        break
        if len(uncoveredRegions[chrom])==0:
            uncoveredRegions.pop(chrom)
    return(uncoveredRegions,amplNames,primersToAmplNames)

def getRegionsUncoveredByDraftExternalPrimers(primersInfo,primersInfoByChrom,outputInternalPrimers,args,refFa):
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
    uncoveredInternalPrimers=deepcopy(outputInternalPrimers)
    regionNameToChrom={}
    amplToStartCoord={}
##    refFa=pysam.FastaFile(args.wholeGenomeRef)
    for chromInt,coords in primersInfoByChrom.items():
        chrom=numToName[chromInt]
        for amplName,parameters in sorted(outputInternalPrimers.items(),
                                          key=lambda item:item[1][3]):
            if chrom!=parameters[3]:
                continue
            elif chrom>parameters[3]:
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
                    if coord[0][0]-100<1:
                        seq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                             chrom,1,coord[1][0]+100)
                    elif coord[1][0]+100>refFa.lengths[refFa.references.index(chromName)]:
                        seq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                             chrom,coord[0][0]-100,refFa.lengths[refFa.references.index(chromName)])
                    else:
                        seq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                             chrom,coord[0][0]-100,coord[1][0]+100)
##                    seq=extractGenomeSeq(refFa,chrom,coord[0][0]-100,coord[1][0]+100)
                    outputExternalPrimers[parameters[2]]=[[leftPrimer,rightPrimer,parameters[2]+'_ext',
                                                          chrom,coord[0][0],coord[1][0],coord[1][0]-coord[0][0]+1,
                                                          coord[0][0]+coord[1][0],coord[1][0]-coord[1][1],
                                                          primersInfo[primerPair][1][0],primersInfo[primerPair][1][1],
                                                          len(leftPrimer),len(rightPrimer),leftGC,rightGC,
                                                          parameters[7]-1-(coord[0][0]+coord[0][1]-1),
                                                          coord[1][0]-coord[1][1]+1-parameters[8]+1,seq]]
                    if amplName in uncoveredInternalPrimers.keys():
                        uncoveredInternalPrimers.pop(amplName)
                    chrTargetSeqStart=parameters[9]-(args.maxExtAmplLen-args.minPrimerShift-args.minPrimerLen)
                    regionNameToChrom[parameters[2]]=chrom
                    amplToStartCoord[parameters[2]]=chrTargetSeqStart
    return(outputExternalPrimers,uncoveredInternalPrimers,regionNameToChrom,amplToStartCoord)

##def writeUnspecificPrimers(primersInfo,rFile,unspecificPrimers):
##    # Write all primers and for primers ]
##    # which give unspecific primers we write primer
##    # with which it form this product
##    wbw=xls.Workbook(rFile)
##    if external:
##        wsw1=wbw.add_worksheet('Draft_Internal_Primers')
##        wsw2=wbw.add_worksheet('Draft_External_Primers')
##    else:
##        wsw2=wbw.add_worksheet('Draft_Internal_Primers')
##        wsw1=wbw.add_worksheet('Draft_External_Primers')
##    wsw1.write_row(0,0,['Primer_Pair','Left_Primer_Start','Left_Primer_Length',
##                        'Right_Primer_End','Right_Primer_Length',
##                        'Left_Primer_Tm','Right_Primer_Tm',
##                        'Amplicon_Length','Primers_Score','Chrom'])
##    rowNum=1
##    for row in oldRows:
##        wsw1.write_row(rowNum,0,row)
##        rowNum+=1
##    wsw2.write_row(0,0,['Primer_Pair','Left_Primer_Start','Left_Primer_Length',
##                       'Right_Primer_End','Right_Primer_Length',
##                       'Left_Primer_Tm','Right_Primer_Tm',
##                       'Amplicon_Length','Primers_Score','Chrom'])
##    rowNum=1
##    for primerPair,info in primersInfo.items():
##        if (goodPrimers==None or
##            goodPrimers and primerPair in goodPrimers):
##            wsw2.write_row(rowNum,0,[primerPair,info[0][0][0],info[0][0][1],
##                                    info[0][1][0],info[0][1][1],
##                                    info[1][0],info[1][1],
##                                    info[2],info[3],info[4]])
##            rowNum+=1
##    wbw.close()

def checkThisPrimerPairForCoveringOtherInputRegions(chromPointRegions,amplBlockStart,amplBlockEnd):
    coveredRegions=[]
    allStarts=[region[1] for region in chromPointRegions.values()]
    allStarts.append(amplBlockStart)
    allStarts.sort()
    for region in sorted(chromPointRegions.values(),
                         key=itemgetter(1))[allStarts.index(amplBlockStart):]:
        chrom,start,end,curRegionName=region
        if start>amplBlockEnd:
            break
        elif start>=amplBlockStart and end<=amplBlockEnd:
            coveredRegions.append(curRegionName)
    return(coveredRegions)

def checkPrimersSpecificity(inputFileBase,primersInfo,
                            wholeGenomeRef,runName,refFa,
                            substNum=1,threads=2,gui=False,
                            maxNonSpecLen=100,
                            maxPrimerNonspec=1000,
                            external=False,varNum=''):
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
    if not os.path.exists(wholeGenomeRef+'.sa'):
        print('WARNING! BWA index is absent for the defined reference genome' '\n' +'Indexing whole genome reference with BWA')
        out=sp.check_output('bwa index '+wholeGenomeRef,shell=True,stderr=sp.STDOUT).decode('utf-8')
    if not os.path.exists(wholeGenomeRef+'.fai'):
        print('WARNING! Samtools index is absent for the defined reference genome' '\n' +'Indexing whole genome reference with samtools')
        out=sp.check_output('samtools faidx '+wholeGenomeRef,shell=True,stderr=sp.STDOUT).decode('utf-8')
    cmd=['bwa','aln','-N','-n',str(args.substNum),
         '-t',str(threads),wholeGenomeRef,seqFile.name]
    with open(bwaResultFileName+'.sai','wb') as file:
        output=sp.call(cmd,stdout=file,stderr=sp.DEVNULL)
    cmd=['bwa','samse','-n','10000000',wholeGenomeRef,
         bwaResultFileName+'.sai',seqFile.name]
    with open(bwaResultFileName+'.sam','wt') as file:
        output=sp.call(cmd,stdout=file,stderr=sp.DEVNULL)
##    out=sp.check_output('bwa aln -N -n '+str(args.substNum)+' -t '+str(threads)+' '+wholeGenomeRef+' '+seqFile.name+' > '+bwaResultFileName+'.sai',shell=True,stderr=sp.STDOUT).decode('utf-8')
##    out=sp.check_output('bwa samse -n 10000000 '+wholeGenomeRef+' '+bwaResultFileName+'.sai'+' '+seqFile.name+' > '+bwaResultFileName+'.sam',shell=True,stderr=sp.STDOUT).decode('utf-8')        
    # Reading BWA output file
    samFile=pysam.AlignmentFile(bwaResultFileName+'.sam')
    # Process SAM-file strings in several threads
    p=ThreadPool(threads)
    print(' Processing SAM-file...')
    logger.info(' Processing SAM-file...')
    totalNumberNonspecificRegions=0
    results=[]
##    refFa=pysam.FastaFile(wholeGenomeRef)
    for read in samFile.fetch():
        results.append(p.apply_async(readBwaFile,(read,maxPrimerNonspec,refFa,
                                                  wholeGenomeRef,primersInfo)))
    wholeWork=len(results)
    doneWork=0
    unspecificPrimers={}
    for res in results:
        # res[0] is a primer name
        res=res.get()
        doneWork+=1
        showPercWork(doneWork,wholeWork,gui)
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
            regions2=primersNonSpecRegions[primerName2]
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
                    if chrom not in regions2.keys():
                        continue
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

def readBwaFile(read,maxPrimerNonspec,refFa,
                wholeGenomeRef,primersInfo=None):
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
        exit(14)
    if primersInfo:
        primersName=read.qname#[:read.qname.rfind('_')]
        for primerPair in primersInfo.keys():
            if primersName in primerPair.split('_'):
                primerNumInPair=primerPair.split('_').index(primersName)
                strand=int((-1)**primerNumInPair)
                pos=int(primersInfo[primerPair][0][primerNumInPair][0])
                try:
                    targetRegion=[primersInfo[primerPair][4],
                                  strand*pos]
                except TypeError:
                    print('ERROR!',primersInfo[primerPair][4],primersInfo[primerPair][0],primerNumInPair)
                    exit(15)
                break
    else:
        targetRegion=[]
    if strand==-1:
        mainMapping=[read.reference_name,strand*(read.pos+len(read.qname))]
    else:
        mainMapping=[read.reference_name,strand*(read.pos+1)]
    if mainMapping!=targetRegion:
        nonSpecRegions=[','.join([read.reference_name,
                                  str(strand*(read.pos+1)),
                                  read.cigarstring,
                                  str(read.get_tag('NM'))])]
    else:
        nonSpecRegions=[]
    if read.has_tag('XA'):
        nonSpecRegions.extend(read.get_tag('XA').split(';')[:-1])
    if len(nonSpecRegions)>maxPrimerNonspec:
        return(read.qname,nonSpecRegions)
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
        if int(pos)<0:
            pos2=-(-int(pos)+regionLen-1)
        else:
            pos2=int(pos)
        if ([chrom,int(pos)]==targetRegion or
            (chrom==targetRegion[0] and
             abs(targetRegion[1]-pos2)<=len(read.seq))):
            continue
        seq=extractGenomeSeq(refFa,wholeGenomeRef,
                             chrom,abs(int(pos)),abs(int(pos))+regionLen-1)
##        attempts=0
##        seq=None
##        while(seq is None):
##            try:
##                seq=refFa.fetch(region=chrom+':'+str(abs(int(pos)))+'-'+str(abs(int(pos))+regionLen-1))
##            except Exception as e:
##                seq=None
##                attempts+=1
##                if attempts>=10:                    
##                    print('ERROR!',e)
##                    logger.error(str(e))
##                    print(refFa.filename)
##                    logger.error(refFa.filename)
##                    exit(16)
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
            try:
                seq=str(Seq.Seq(seq).reverse_complement()).upper()
            except ValueError:
                print('ERROR! Unknown DNA symbols in the sequence extracted from the reference genome!')
                print(seq)
                print(chrom,pos,regionLen)
                exit(49)
        else:
            seq=seq.upper()
        # If there is some insertions or deletion in found region
        if 'I' in cigar or 'D' in cigar:
            # We need to align its sequence with primer sequence
            align=pairwise2.align.globalxx(read.seq,seq)
            # Check that 3'-ends of primers are identical
            ## and one of two nucleotides before 3'-ends are identical, too
            if (align[0][0][-1]==align[0][1][-1] and
                (align[0][0][-2]==align[0][1][-2] or
                 align[0][0][-3]==align[0][1][-3])):
                # Then we consider this region as a non-specific for this primer
                primerNonSpecRegions.append([chrom,
                                             regStrand,
                                             abs(int(pos)),
                                             regionLen])
        else:
            if seq[-1]==read.seq[-1] and (seq[-2]==read.seq[-2] or
                                         seq[-3]==read.seq[-3]):
                # Then we consider this region as a non-specific for this primer
                primerNonSpecRegions.append([chrom,
                                             regStrand,
                                             abs(int(pos)),
                                             regionLen])
    return(read.qname,primerNonSpecRegions)

def getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,maxNonSpecLen=100,
                                            threads=2,gui=False):
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
        showPercWork(i+1,allWork,gui)
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

def removeBadPrimerPairs(primersInfoByChrom,primersInfo,goodPrimers,primersToAmplNames,amplNames):
    newPrimersInfoByChrom={}
    for chromInt,primers in primersInfoByChrom.items():
        chrom=numToName[chromInt]
        for primerPairName,primers in primers.items():
            if primerPairName in goodPrimers:
                if chromInt not in newPrimersInfoByChrom.keys():
                    newPrimersInfoByChrom[chromInt]={primerPairName:primers}
                elif primerPairName not in newPrimersInfoByChrom[chromInt].keys():
                    newPrimersInfoByChrom[chromInt][primerPairName]=primers
                else:
                    print('ERROR! Pair of primers is repeated in the primersInfoByChrom!')
                    print(chrom,primerPairName)
                    logger.error('Pair of primers is repeated in the primersInfoByChrom!')
                    logger.error(chrom,primerPairName)
                    exit(17)
            else:
                # If draft file contained some primer pairs
                # that do not overlap any target region,
                # we can simply skip them
                if primerPairName in primersToAmplNames.keys():
                    amplNamesToDelete=primersToAmplNames[primerPairName]
                    for amplName in amplNamesToDelete:
                        curRegionName=amplName[:amplName.rfind('_')]
                        amplNames[curRegionName].remove(amplName)
                primersInfo.pop(primerPairName)
    primersInfoByChrom=deepcopy(newPrimersInfoByChrom)
    return(primersInfoByChrom,primersInfo,amplNames)

def checkThatAllInputRegionsCovered(amplNames,allRegions,
                                    regionNameToChrom,regionsCoords,
                                    filterMessage='by specificity',skipUndesigned=False):
    someInputRegionUncovered=False
    newAmplNames={}
    for curRegionName,ampls in sorted(amplNames.items(),
                                      key=lambda item:splitNameForSorting(item[0])):
        if len(ampls)==0:
            if skipUndesigned:
                try:
                    print('WARNING! For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                except KeyError:
                    print('KeyError:',curRegionName)
                    print(curRegionName in regionNameToChrom.keys())
                    print(regionNameToChrom[curRegionName] in allRegions.keys())
                    print(curRegionName in allRegions[regionNameToChrom[curRegionName]].keys())
                    exit(18)
                logger.warn('For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                regionsCoords[nameToNum[regionNameToChrom[curRegionName]]].remove(allRegions[regionNameToChrom[curRegionName]][curRegionName][1])
                if len(regionsCoords[nameToNum[regionNameToChrom[curRegionName]]])==0:
                    regionsCoords.pop(nameToNum[regionNameToChrom[curRegionName]])
                allRegions[regionNameToChrom[curRegionName]].pop(curRegionName)
                if len(allRegions[regionNameToChrom[curRegionName]])==0:
                    allRegions.pop(regionNameToChrom[curRegionName])
                regionNameToChrom.pop(curRegionName)
            else:
                newAmplNames[curRegionName]=ampls
                print('ERROR! For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                logger.error('For input region '+curRegionName+' ('+' '.join(list(map(str,allRegions[regionNameToChrom[curRegionName]][curRegionName])))+') no primers left after filtering primers '+filterMessage+'! Try to increase -primernum1 parameter. Or if you have already tried, use less stringent parameters.')
                someInputRegionUncovered=True
        else:
            newAmplNames[curRegionName]=ampls
    if someInputRegionUncovered:
        exit(19)
    return(newAmplNames,allRegions,regionNameToChrom,regionsCoords)

def analyzePrimersForCrossingSNP(primersInfo,threads,
                                 dbSnpVcfFile,
                                 end3Len=None,gui=False):
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
            exit(20)
        strand1=1; strand2=-1
        if primer1!='':
            results.append(p.apply_async(checkPrimerForCrossingSNP,(primer1,chrom,start1,end1,strand1,
                                                                    dbSnpVcfFile,args.snpFreq,
                                                                    end3Len)))
        if primer2!='':
            results.append(p.apply_async(checkPrimerForCrossingSNP,(primer2,chrom,start2,end2,strand2,
                                                                    dbSnpVcfFile,args.snpFreq,
                                                                    end3Len)))
    primersCoveringSNPs=[]
    wholeWork=len(results)
    done=0
    for res in results:
        primer,result=res.get()
        if result:
            primersCoveringSNPs.append(primer)
        done+=1
        showPercWork(done,wholeWork,gui)
    p.close()
    p.join()
    primerPairsNonCoveringSNPs=[]
    for primers,info in primersInfo.items():
        primer1,primer2=primers.split('_')
        if primer1 not in primersCoveringSNPs and primer2 not in primersCoveringSNPs:
            primerPairsNonCoveringSNPs.append(primers)
    print("\n # Number of primers that do not overlap with high-frequent SNPs: "+str(len(primerPairsNonCoveringSNPs)))
    logger.info("\n # Number of primers that do not overlap with high-frequent SNPs: "+str(len(primerPairsNonCoveringSNPs)))
    return(primerPairsNonCoveringSNPs,primersCoveringSNPs)    

# checkPrimerForCrossingSNP checks, if primer crosses some SNPs with high frequence in population
def checkPrimerForCrossingSNP(primer,chrom,start,end,strand,
                              dbSnpVcfFile,freq=0.1,end3Len=None):
    vcf=pysam.VariantFile(dbSnpVcfFile)
    # end3Len sets number of nucleotides from the 3'-end that we will check
    if end3Len!=None:
        if strand>0:
            end3_start=max(start,end-end3Len+1)
            end3=list(range(end3_start,end+1))
        else:
            end3_end=min(end,start+end3Len-1)
            end3=list(range(start,end3_end+1))
    else:
        end3=list(range(start,end+1))
    # Previously, we only filtered out primers that crossed SNP by 3'-end
    # Now, we filter primers that cross SNP by any of its sequence
    hfSnpFound=False
    snps=[]
    snpsEnd3=[]
    try:
        for snp in vcf.fetch(chrom,start-1,end):
            if float(snp.info['CAF'][0])<1-freq:
                snps.append(snp)
                if snp.pos in end3:
                    snpsEnd3.append(snp)
    except ValueError:
        try:
            for snp in vcf.fetch(chrom.replace('chr',''),start-1,end):
                if float(snp.info['CAF'][0])<1-freq:
                    snps.append(snp)
                    if snp.pos in end3:
                        snpsEnd3.append(snp)
        except:
            print("ERROR (49)! Unknown error with getting SNPs from the following file. Possibly you haven't made tbi-file with tabix:")
            logger.error("(49)! Unknown error with getting SNPs from the following file. Possibly you haven't made tbi-file with tabix:")
            print('File:',dbSnpVcfFile)
            logger.error('File: '+dbSnpVcfFile)
            print('Chromosome:',chrom)
            logger.error('Chromosome: '+chrom)
            print('Start:',start-1)
            logger.error('Start: '+str(start-1))
            print('End:',end)
            logger.error('End: '+str(end))
            exit(49)
    if len(snps)>1 or len(snpsEnd3)>0:
        hfSnpFound=True
    return(primer,hfSnpFound)

# joinAmpliconsToAmplifiedBlocks joins neighbourhing or overlapping amplicons to blocks
# Minimal path is a path of one primer:
## [coord1,primer_name,coord2]
# chrom here is a chromosome number
def joinAmpliconsToBlocks(chromRegionsCoords,chromPrimersInfoByChrom,
                          maxAmplLen=100,chromInt=None,
                          returnVariantsNum=10):
    chrom=numToName[chromInt]
    # First, split all regions on this chromosome onto blocks, elements of which cannot be joined into one amplificated block
    blocks=[[chromRegionsCoords[0]]]
    coordToBlock={} # converts coordinate into the number of block
    for coord in chromRegionsCoords[1:]:
        if coord>=blocks[-1][-1]+maxAmplLen*2:
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
    comparisonNum=0
    for i,(primerPairName1,primers1) in enumerate(sorted(chromPrimersInfoByChrom.items(),
                                                         key=lambda item:(item[1][0][0],
                                                                          item[1][1][0]))):
        primerName1,primerName2=primerPairName1.split('_')
        amplBlockStart1=primers1[0][0]+primers1[0][1]
        amplBlockEnd1=primers1[1][0]-primers1[1][1]
        # Get block number
        coords=deepcopy(chromRegionsCoords)
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
        if amplBlockEnd1 in blocks[blockNum]:
            lastMutNum1=sorted(blocks[blockNum]).index(amplBlockEnd1)
        else:
            coords=deepcopy(blocks[blockNum])
            coords.append(amplBlockEnd1)
            lastMutNum1=sorted(coords).index(amplBlockEnd1)-1
        if amplBlockStart1 in blocks[blockNum]:
            firstMutNum1=sorted(blocks[blockNum]).index(amplBlockStart1)
        else:
            coords=deepcopy(blocks[blockNum])
            coords.append(amplBlockStart1)
            firstMutNum1=sorted(coords).index(amplBlockStart1)
        nextMutNotCovered=True
        prevMutNotCovered=True
        for primerPairName2,primers2 in sorted(chromPrimersInfoByChrom.items(),
                                               key=lambda item:(item[1][0][0],
                                                                item[1][1][0]))[i+1:]:
            if primerPairName1==primerPairName2:
                continue
            comparisonNum+=1
            amplBlockStart2=primers2[0][0]+primers2[0][1]
            amplBlockEnd2=primers2[1][0]-primers2[1][1]
            nextMut=blocks[blockNum][min(lastMutNum1+1,len(blocks[blockNum])-1)]
            prevMut=blocks[blockNum][max(firstMutNum1-1,0)]
            if (not lastMut and
                amplBlockStart2<=nextMut<=amplBlockEnd2 and
                amplBlockEnd1<amplBlockStart2+args.maxoverlap):
                # Calculate weight of this edge - distance between amplicons
                weight=maxAmplLen-min(50,primers2[0][0]-primers1[1][0]) # if distance is too large (>50), leave it as 50
                blockGraph.add_edge(primerPairName1,primerPairName2,attr_dict={'weight':weight})
                nextMutNotCovered=False
            elif (not firstMut and
                  amplBlockStart2<=prevMut<=amplBlockEnd2 and
                  amplBlockEnd2<amplBlockStart1+args.maxoverlap):
                # Calculate weight of this edge - distance between amplicons
                weight=maxAmplLen-min(50,primers1[0][0]-primers2[1][0]) # if distance is too large (>50), leave it as 50
                blockGraph.add_edge(primerPairName1,primerPairName2,attr_dict={'weight':weight})
                prevMutNotCovered=False
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
                exit(24)
            finalShortestPaths=[]
            for path in paths:
                finalShortestPaths.append(path[1:-1])
        else:
            try:
                path=tuple(nx.algorithms.shortest_paths.generic.shortest_path(blockGraph,firstNodes[i],lastNodes[i]))
            except nx.exception.NetworkXNoPath as e:
                print('ERROR! Too low value of maximal overlap (-maxoverlap) or of initially designed primers (-primernum): '+str(args.maxoverlap)+' and '+str(args.primernum1)+'. Try to increase one of them')
                logger.error(' Too low value of maximal overlap (-maxoverlap) or of initially designed primers (-primernum): '+str(args.maxoverlap)+' and '+str(args.primernum1)+'. Try to increase one of them')
                print(str(lastNodes[i])+' is not reachable from '+str(firstNodes[i]))
                logger.error(str(lastNodes[i])+' is not reachable from '+str(firstNodes[i]))
                logger.error(str(blockGraph.edges()))
                logger.error(str(chromPrimersInfoByChrom))
                exit(25)
            g1=blockGraph.subgraph(path)
            shortestPaths={path:len(path)}
            minPathLen=len(path)
            maxAnalysisVars=returnVariantsNum*2
            analysisNum=0
            weightOff=False
            for path in nx.shortest_simple_paths(blockGraph,firstNodes[i],lastNodes[i]):
                analysisNum+=1
                shortestPaths[tuple(path)]=len(path)
                if len(shortestPaths.keys())>=returnVariantsNum*2:
                    break
                if analysisNum>maxAnalysisVars:
                    break
            finalShortestPaths=[]
            j=0
            for path,value in sorted(shortestPaths.items(),key=itemgetter(1)):
                if 'end' in path[1:-1]:
                    continue
                finalShortestPaths.append(path[1:-1])
                j+=1
                if j>=args.returnVariantsNum:
                    break
        blocksFinalShortestPaths.append(finalShortestPaths)
    return(chromInt,blocksFinalShortestPaths)

# allRegionsAmplifiedBlocks contains chromosome numbers
def getBestPrimerCombinations(allRegionsAmplifiedBlocks,
                              primersInfo,
                              unspecificPrimers,
                              returnVarNum=1,
                              triesToGetCombination=1000,
                              totalMultiplexVariants=1000,
                              gui=False):
    print(' Searching for interactions between the amplified blocks...')
    logger.info(' Searching for interactions between the amplified blocks...')
    # Save interaction numbers for all possible pairs of blocks
    allBlockPairValues={}
    # Total number of blocks
    totalBlockNum=0
    # Maximal number of block variants
    maxBlockVarNum=0
    # Go through all blocks
    # Go through all chromosomes
    showPercWork(0,len(allRegionsAmplifiedBlocks.keys()),gui)
    for i,chromInt1 in enumerate(sorted(allRegionsAmplifiedBlocks.keys())):
        totalBlockNum+=len(allRegionsAmplifiedBlocks[chromInt1])
        # Go through all blocks of current chromosome
        for k,block1 in enumerate(allRegionsAmplifiedBlocks[chromInt1]):
            if len(block1)>maxBlockVarNum:
                maxBlockVarNum=len(block1)
            # Go through all chromosomes
            for j,chromInt2 in enumerate(sorted(allRegionsAmplifiedBlocks.keys())):
                # If chrom2 is sorted earlier than chrom1, skip chrom2
                if j<i:
                    continue
                # If it is the last block of chrom1, we skip this chromosome
                # Remove variant of the last block of chrom1
                # Because we changed functions of calculating interNum
                # and now we need to compare each block with itself
                # If it is the same chromosome
                else:
                    if chromInt1==chromInt2:
                        startNum=k
                    else:
                        startNum=0
                    # Go through all blocks of the 2nd blocks to compare
                    for z,block2 in enumerate(allRegionsAmplifiedBlocks[chromInt2][startNum:]):
                        if len(block2)>maxBlockVarNum:
                            maxBlockVarNum=len(block2)
                        # Go through all variants of current block1
                        for q,block1_var in enumerate(block1):
                            # Go through all variants of current block2
                            for w,block2_var in enumerate(block2):
                                # Get number of interacting primers
                                interNum=comparePrimersOfTwoBlocks(block1_var,block2_var,unspecificPrimers)
                                block_ID1='_'.join([str(chromInt1),str(k),str(q)])
                                block_ID2='_'.join([str(chromInt2),str(startNum+z),str(w)])
                                if chromInt1 not in allBlockPairValues.keys():
                                    allBlockPairValues[chromInt1]={}
                                if k not in allBlockPairValues[chromInt1].keys():
                                    allBlockPairValues[chromInt1][k]={}
                                if q not in allBlockPairValues[chromInt1][k].keys():
                                    allBlockPairValues[chromInt1][k][q]={}
                                allBlockPairValues[chromInt1][k][q][block_ID2]=interNum
                                if chromInt2 not in allBlockPairValues.keys():
                                    allBlockPairValues[chromInt2]={}
                                if startNum+z not in allBlockPairValues[chromInt2].keys():
                                    allBlockPairValues[chromInt2][startNum+z]={}
                                if w not in allBlockPairValues[chromInt2][startNum+z].keys():
                                    allBlockPairValues[chromInt2][startNum+z][w]={}
                                allBlockPairValues[chromInt2][startNum+z][w][block_ID1]=interNum
        showPercWork(i+1,len(allRegionsAmplifiedBlocks.keys()),gui)
    print('\n Searching for the best primer pair combinations...')
    logger.info(' Searching for the best primer pair combinations...')
    # Create matrix of sorted values for each of block variants:
    # [[block1-best,block2-best,block3-best,...blockN-best]
    #  [block1-next,block2-next,block3-next,...blockN-next]
    #   ...
    #  [block1-next,block2-next,block3-next,...blockN-next]]
    blockVarValues=numpy.full((maxBlockVarNum,
                               totalBlockNum),None)
    # Create matrix of block variant names sorted by values:
    # [[block1-best,block2-best,block3-best,...blockN-best]
    #  [block1-next,block2-next,block3-next,...blockN-next]
    #   ...
    #  [block1-next,block2-next,block3-next,...blockN-next]]
    blockVarNames=numpy.full((maxBlockVarNum,
                              totalBlockNum),None)
    # Stores positional number of current block among all blocks
    currentBlockNum=0
    # Go through all chromosomes
    for chromInt in sorted(allRegionsAmplifiedBlocks.keys()):
        # Go through all blocks of current chromosome
        for k,block in enumerate(allRegionsAmplifiedBlocks[chromInt]):
            for j,blockVarNum in enumerate(sorted(allBlockPairValues[chromInt][k],
                                                  key=lambda blockVarNum:sum(allBlockPairValues[chromInt][k][blockVarNum].values()))):
                # Save value for current block variant
                blockVarValues[j,currentBlockNum]=sum(allBlockPairValues[chromInt][k][blockVarNum].values())
                # Create block variant name
                block_ID='_'.join([str(chromInt),str(k),str(blockVarNum)])
                # Save name of current block into the matrix
                blockVarNames[j,currentBlockNum]=block_ID
            currentBlockNum+=1
    # Stores all best combinations
    combinations={}
    # Calculate score for the best variant
    score=sum(blockVarValues[0,])
    # Create first best combination
    combinations[tuple([0]*totalBlockNum)]=score
    # Stores number of created combinations
    totalCombNum=0
    # Create combinations while we haven't created enough number
    # Or while it wasn't interrupted in cycle
    # At the 1st step we create two time more combinations
    while (totalCombNum<min(triesToGetCombination,totalMultiplexVariants)):
        # Stores all possible next combinations with sum values
        nextPossibleCombinations={}
        # Go through all created combinations and find block
        # with minimal difference between current and the next variant
        for combination,score in combinations.items():
            # Go through all blocks
            for i,blockVarNumInSorted in enumerate(combination):
                # If current block variant number is not None
                # i.e. it exists
                if (blockVarValues.shape[0]>=blockVarNumInSorted+2 and
                    blockVarValues[blockVarNumInSorted+1,i]!=None):
                    # Change only one block variant number
                    newCombination=combination[:i]+tuple([combination[i]+1])+combination[i+1:]
                    if tuple(newCombination) not in combinations.keys():
                        diff=blockVarValues[blockVarNumInSorted+1,i]-blockVarValues[blockVarNumInSorted,i]
                        nextPossibleCombinations[tuple(newCombination)]=score+diff
        # If there is no next possible combinations, interrupt cycles
        if len(nextPossibleCombinations)==0:
            break
        # Go through all saved possible combinations and take with minimal score
        minScore=min(nextPossibleCombinations.values())
        for combination,value in sorted(nextPossibleCombinations.items(),
                                        key=lambda item:item[1]):
            if value>minScore:
                break
            # Save combinations with minimal scores
            combinations[combination]=value
        totalCombNum=len(combinations.keys())
        showPercWork(totalCombNum,min(triesToGetCombination,totalMultiplexVariants),gui)
    print()
    print(' # Total number of combinations selected for the subsequent deeper analysis:',len(combinations))
    logger.info(' # Total number of combinations selected for the subsequent deeper analysis: '+str(len(combinations)))
    print(' # Minimal score of the combinations before recalculation:',min(combinations.values()))
    logger.info(' # Minimal score of the combinations analyzed before recalculation: '+str(min(combinations.values())))
    print(' # Maximal score of the combinations analyzed before recalculation:',max(combinations.values()))
    logger.info(' # Maximal score of the combinations analyzed before recalculation: '+str(max(combinations.values())))
    print(' Recalculating scores for the selected combinations and taking best ones...')
    logger.info(' Recalculating scores for the selected combinations and taking best ones...')
    # At the 2nd step we take only best combinations
    # but by the more complicated comparison
    # Previously we considered interactions with primers from all variants
    # Now we will consider only interactions with primers from current combination
    # This will give more precise evaluation of combination score
    # Go through all combinations and their scores
    for combNum,(combination,score) in enumerate(combinations.items()):
        # combination is a tuple of block variant numbers in sorted list
        # We need to count interactions between all primers of these block variants
        # Stores list of all blocks of current combination in the format
        # that is necessary for comparePrimersOfTwoBlocks()
        combAllBlocks=[]
        for i,blockVarNameNum in enumerate(combination):
            chromInt,blockNum,blockVarNum=blockVarNames[blockVarNameNum,i].split('_')
            blockPrimers=allRegionsAmplifiedBlocks[int(chromInt)][int(blockNum)][int(blockVarNum)]
            combAllBlocks.append(blockPrimers)
        newScore=0
        # Count interactions between all primers of all blocks
        for i,block1 in enumerate(combAllBlocks):
            for block2 in combAllBlocks[i+1:]:
                newScore+=comparePrimersOfTwoBlocks(block1,block2,unspecificPrimers)
        combinations[combination]=newScore
        showPercWork(combNum+1,len(combinations),gui)
    print('\n # Minimal score of the combinations analyzed:',min(combinations.values()))
    logger.info(' # Minimal score of the combinations analyzed: '+str(min(combinations.values())))
    print(' # Maximal score of the combinations analyzed:',max(combinations.values()))
    logger.info(' # Maximal score of the combinations analyzed: '+str(max(combinations.values())))
    # Select necessary number of combinations
    # Contains combinations in the necessary format for the subsequent steps
    outCombinations=[]
    for i,combination in enumerate(sorted(combinations.keys(),
                                          key=lambda key:combinations[key])):
        if i+1<=returnVarNum:
            outCombinations.append({})
        else:
            break
        primerPairNum=0
        for j,blockVarNameNum in enumerate(combination):
            chromIntStr,blockNum,blockVarNum=blockVarNames[blockVarNameNum,j].split('_')
            blockPrimers=allRegionsAmplifiedBlocks[int(chromIntStr)][int(blockNum)][int(blockVarNum)]
            if chromIntStr not in outCombinations[-1].keys():
                outCombinations[-1][chromIntStr]={}
            for primerPair in blockPrimers:
                try:
                    outCombinations[-1][chromIntStr][primersInfo[primerPair][0][0][0]]=primerPair.split('_')
                except KeyError:
                    print('ERROR (59)! Unknown error with primerPair:')
                    logger.error('(59)! Unknown error with primerPair:')
                    print('PrimerPair: '+str(primerPair))
                    logger.error('PrimerPair: '+str(primerPair))
                    print(blockPrimers)
                    logger.error(str(blockPrimers))
                    print(allRegionsAmplifiedBlocks)
                    exit(59)
                primerPairNum+=1
        print(' # Number of primer pairs and score for the multiplex variant',i+1,'-',primerPairNum,'('+str(combinations[combination])+')')
        logger.info(' # Number of primer pairs and score for the multiplex variant '+str(i+1)+': '+str(primerPairNum)+' ('+str(combinations[combination])+')')
    return(outCombinations)

def getEdgesSumOfWeight(g,nodes):
    sumWeight=0
    for i,node in enumerate(nodes[1:-1]):
        sumWeight+=g[node][nodes[i]]['weight']
    return(sumWeight)

# Returns number of interactions for two blocks
def comparePrimersOfTwoBlocks(block1,block2,unspecificPrimers):
    interactions=0
    for primerPair1 in block1:
        for primerPair2 in block2:
            for primer1 in primerPair1.split('_'):
                for primer2 in primerPair2.split('_'):
                    if ('_'.join([primer1,primer2]) in unspecificPrimers or
                        '_'.join([primer2,primer1]) in unspecificPrimers):
                        interactions+=1
                    # Check only for 6 nucleotides that can give hybridized 3'-end
                    elif revComplement(primer1[-6:]) in primer2:
                        interactions+=1
                    elif revComplement(primer2[-6:]) in primer1:
                        interactions+=1
    return(interactions)

def makeFinalMultiplexes(initialGraph,multiplexes=[],multNum=2,functionFirstCall=False):
    # We try to make cliques from it
    if len(initialGraph)==0:
        return(multiplexes)
    elif len(initialGraph.nodes())==1:
        return(multiplexes+[initialGraph.nodes()])
    elif (len(initialGraph.nodes())==multNum and
          len(multiplexes)==0):
        return([[x] for x in initialGraph.nodes()])
    graph=deepcopy(initialGraph)
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
        exit(26)
    if len(cls[0])>maxNodesNum:
        # Sort by number of nodes to which this node is linked by edges
        cls[0]=sorted(cls[0],key=graph.degree)
        cls[0]=cls[0][:-(len(cls[0])-maxNodesNum)]
    multiplexes.append(cls[0])
    graph.remove_nodes_from(cls[0])
    if len(graph)>0 and multNum-1>0:
        multiplexes=makeFinalMultiplexes(graph,multiplexes,multNum-1,False)
    if functionFirstCall:
        for multiplex in multiplexes:
            for multAmplicon in multiplex:
                if multAmplicon in graph.nodes():
                    graph.remove_node(multAmplicon)
    # If there are some nodes that were not sorted to any of multiplexes
    # This is the first call of this function
    if len(graph.nodes())>0 and functionFirstCall:
        currentMultNumToAdd=0
        # Go through all undorted nodes
        for leftUnsortedNode in graph.nodes():
            multNumsToWhichWeTried=[]
            neighbours=initialGraph.neighbors(leftUnsortedNode)
            # Go through all multiplexes and try to add the unsorted node to this multiplex
            while(len(multNumsToWhichWeTried)<multNum):
                thisAmpliconFitsThisMultiplex=True
                # If current unsorted amplicon is linked to all amplicons of current multiplex,
                # add this amplicons to this multiplex
                for multAmplicon in multiplexes[currentMultNumToAdd]:
                    if multAmplicon not in neighbours:
                        thisAmpliconFitsThisMultiplex=False
                        break
                if thisAmpliconFitsThisMultiplex:
                    multiplexes[currentMultNumToAdd].append(leftUnsortedNode)
                    break
                multNumsToWhichWeTried.append(currentMultNumToAdd)
                currentMultNumToAdd+=1
                if currentMultNumToAdd>multNum-1:
                    currentMultNumToAdd=0
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
        cls=makeFinalMultiplexes(localMultiplexNums,[],len(mults),True)
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
            leftNodesGraph=deepcopy(localMultiplexNums)
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
            exit(27)
    return(multiplexes)

# Function that checks that new primer fits one of the multiplexes
## primers is a list of two primers and their amplicon coordinates [primer1, primer2, amplStart, amplEnd]
## nums are numbers of multiplexes to which we try to put new primer
## globalMultiplexNums are lists of multiplexes with primers and their amplicon coordinates [primer1, primer2, amplStart, amplEnd]
## unspecificPrimers is a list of primer pairs that form unspecific products 
def checkPrimersFit(primers,primersToCompare,
                    minmultdimerdg1=-6,minmultdimerdg2=-10,
                    maxIntersectionOfPrimers=5,
                    unspecificPrimers=[],
                    mv_conc=50,dv_conc=3,
                    dntp_conc=0.8,dna_conc=250,
                    leftAdapter=None,rightAdapter=None):
    leftPrimer,rightPrimer=primersToCompare[0:2]
    chrom,amplStart,amplEnd=primersToCompare[3:6]
    # We need to check the following:
    ## Current primers pair does not overlap with primersToCompare (the most important!)
    ## Current primers pair does not form any secondary structure (important, if with 5'-overhang)
    ## Current primers pair does not form any unspecific product
    ## Maybe, GC-content difference
    # Overlapping is the most important point of checking
    # Here we compare chromosome names
    if chrom==primers[2]:
        inter=set(range(primers[3],primers[4])).intersection(list(range(amplStart,amplEnd)))
        if len(inter)>maxIntersectionOfPrimers:
            # If there is intersection of length more than 5 bp, these primers pair do not correspond each other
            return(False,'Amplicons overlap')
    # Secondary structure with 5'-overhang
    maxdG=minmultdimerdg1
    maxdG2=minmultdimerdg2
    if leftAdapter:
        leftAdapter=leftAdapter.upper()
        leftPrimerToCheck1=leftAdapter+primers[0]
        leftPrimerToCheck2=leftAdapter+leftPrimer
    else:
        leftPrimerToCheck1=primers[0]
        leftPrimerToCheck2=leftPrimer
    if rightAdapter:
        rightAdapter=rightAdapter.upper()
        rightPrimerToCheck1=rightAdapter+primers[1]
        rightPrimerToCheck2=rightAdapter+rightPrimer
    else:
        rightPrimerToCheck1=primers[1]
        rightPrimerToCheck2=rightPrimer
    dG1=primer3.calcHeterodimer(leftPrimerToCheck1,leftPrimerToCheck2,
                                mv_conc=mv_conc,dv_conc=dv_conc,
                                dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1000
    dG2=primer3.calcHeterodimer(leftPrimerToCheck1,rightPrimerToCheck2,
                                mv_conc=mv_conc,dv_conc=dv_conc,
                                dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1000
    dG3=primer3.calcHeterodimer(rightPrimerToCheck1,leftPrimerToCheck2,
                                mv_conc=mv_conc,dv_conc=dv_conc,
                                dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1000
    dG4=primer3.calcHeterodimer(rightPrimerToCheck1,rightPrimerToCheck2,
                                mv_conc=mv_conc,dv_conc=dv_conc,
                                dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1000
    try:
        if (calcThreeStrikeEndDimer(leftPrimerToCheck1,leftPrimerToCheck2,
                                    mv_conc,dv_conc,
                                    dntp_conc,dna_conc)<maxdG
            or dG1<maxdG2):
            return(False,'Heterodimer of F-primer with F-primer')
    except TypeError as e:
        print('ERROR! Unknown error:',e)
        logger.error('Unknown error: '+str(e))
        print(calcThreeStrikeEndDimer(leftPrimerToCheck1,leftPrimerToCheck2,mv_conc,dv_conc,dntp_conc,dna_conc))
        print(maxdG)
        logger.error(calcThreeStrikeEndDimer(leftPrimerToCheck1,leftPrimerToCheck2,mv_conc,dv_conc,dntp_conc,dna_conc))
        logger.error(maxdG)
        exit(51)
    if (calcThreeStrikeEndDimer(leftPrimerToCheck1,rightPrimerToCheck2,
                                mv_conc,dv_conc,
                                dntp_conc,dna_conc)<maxdG
        or dG2<maxdG2):
        return(False,'Heterodimer of F-primer with R-primer')
    if (calcThreeStrikeEndDimer(rightPrimerToCheck1,leftPrimerToCheck2,
                                mv_conc,dv_conc,
                                dntp_conc,dna_conc)<maxdG
        or dG3<maxdG2):
        return(False,'Heterodimer of R-primer with F-primer')
    if (calcThreeStrikeEndDimer(rightPrimerToCheck1,rightPrimerToCheck2,
                                mv_conc,dv_conc,
                                dntp_conc,dna_conc)<maxdG
        or dG4<maxdG2):
        return(False,'Heterodimer of R-primer with R-primer')
    # Unspecific products
    if len(unspecificPrimers)>0:
        if ('_'.join([primers[0],leftPrimer]) in unspecificPrimers
            or '_'.join([primers[0],rightPrimer]) in unspecificPrimers
            or '_'.join([primers[1],leftPrimer]) in unspecificPrimers
            or '_'.join([primers[1],rightPrimer]) in unspecificPrimers):
            return(False,'Unspecific product')
    return(True,None)

# Function that checks if two primers can hybridize with hybridized 3'-end
# Start from 4 nucleotides
def calcThreeStrikeEndDimer(primer1,primer2,
                            mv_conc=50,dv_conc=3,
                            dntp_conc=0.8,dna_conc=250,
                            startEndLen=4):
    # Get the longest region of primer1 3'-end that can hybridize to primer2
    endLen=startEndLen
    if revComplement(primer1[-endLen:]) in primer2:
        while(endLen<len(primer1)):
            endLen+=1
            if revComplement(primer1[-endLen:]) not in primer2[1:]:
                endLen-=1
                break
        # Get dG for the longest region found
        # We take -1 and +1 nucleotides from the hybridized region
        thermoStart1=min(len(primer1),endLen+1)
        thermoStart2=max(0,primer2.index(revComplement(primer1[-endLen:]))-1)
        thermoEnd2=min(len(primer2),primer2.index(revComplement(primer1[-endLen:]))+endLen+1)
        dG1=primer3.calcHeterodimer(primer1[-thermoStart1:],primer2[thermoStart2:thermoEnd2],
                                    mv_conc=mv_conc,dv_conc=dv_conc,
                                    dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1000
    else:
        dG1=0
    if primer1!=primer2:
        # Get the longest region of primer2 3'-end that can hybridize to primer1
        endLen=startEndLen
        if revComplement(primer2[-endLen:]) in primer1[1:]:
            while(endLen<len(primer2)):
                endLen+=1
                if revComplement(primer2[-endLen:]) not in primer1:
                    endLen-=1
                    break
            # Get dG for the longest region found
            # We take -1 and +1 nucleotides from the hybridized region
            thermoStart1=min(len(primer2),endLen+1)
            thermoStart2=max(0,primer1.index(revComplement(primer2[-endLen:]))-1)
            thermoEnd2=min(len(primer1),primer1.index(revComplement(primer2[-endLen:]))+endLen+1)
            dG2=primer3.calcHeterodimer(primer2[-thermoStart1:],primer1[thermoStart2:thermoEnd2],
                                        mv_conc=mv_conc,dv_conc=dv_conc,
                                        dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1000
        else:
            dG2=0
    else:
        dG2=0
    return(min(dG1,dG2))

# Function that checks if primer can form hairpin with hybridized 3'-end
# Start from 3 nucleotides
def calcThreeStrikeEndHairpin(seq,
                              nucs=3,minLoop=3,
                              mv_conc=50,dv_conc=3,
                              dntp_conc=0.8,dna_conc=250):
    seq=seq.upper()
    # Create reverse complement sequence for the sequence
    seq_rc=revComplement(seq)
    # Take last 3 nucleotides of 3'-end of sequence
    end3=seq[-nucs:]
    # Check if end3 has reverse complement sequence
    # in other part of sequence with distance minimum 3 nucs
    if end3 in seq_rc[len(end3)+minLoop:]:
        # Then trying to extend this end3
        # Var for storing extension length
        extension=0
        while(True):
            # Increase extension
            extension+=1
            if nucs+extension>len(seq)-(nucs+minLoop+extension):
                end3=seq[-(nucs+extension):]
                break
            # Extend end3
            end3=seq[-(nucs+extension):]
            if end3 not in seq_rc[len(end3)+minLoop:]:
                # We do not remove 1st nucleotide because later 
                # it will be necessary for dimer dG estimation
                break
        # Extract rev-compl sequence for it
        # First, get start of it in seq_rc
        start=seq_rc[len(end3)-1+minLoop:].index(end3[1:])+len(end3)+minLoop-1
        # Get end
        end=start+len(end3)
        # Extract from seq subsequence that will hybridize with end3
        subSeq=seq[-end:-(start-1)]
        # When we have the longest part of end3 that is revComplement
        # to some part of seq,
        # approximately estimate its strength by calculating dG for dimer
##        print(subSeq)
##        print(end3)
        dG=primer3.calcHeterodimer(subSeq,end3,
                                   mv_conc=mv_conc,dv_conc=dv_conc,
                                   dntp_conc=dntp_conc,dna_conc=dna_conc).dg/1500
        # We divide it on 2000 instead of 1000, because it's only approximation
        return(dG)
    else:
        return(0)

def extractGenomeSeq(refFa,wholeGenomeRef,
                     chromName,start,end):
    attempts=0
    seq=None
    while(seq is None):
        try:
            seq=refFa.fetch(region=chromName+':'+str(start)+'-'+str(end))
        except:
            seq=None
            attempts+=1
            # If there are some problems with getting regions seq,
            # re-open whole genome reference file
            refFa=pysam.FastaFile(wholeGenomeRef)
            if attempts>=10:                    
                print('ERROR: could not extract genome sequence!',)
                logger.error('Could not extract genome sequence!')
                print(refFa.filename)
                logger.error(refFa.filename)
                print(chromName,start,end)
                logger.error(chromName+':'+str(start)+'-'+str(end))
                print('Chromosome name to number converter:',nameToNum)
                logger.error('Chromosome name to number converter: '+str(nameToNum.items()))
                print('Chromosome number to name converter:',numToName)
                logger.error('Chromosome number to name converter: '+str(numToName.items()))
                exit(28)
    return(seq.upper())

def showPercWork(done,allWork,gui=False):
    percDoneWork=round((done/allWork)*100,2)
    if (percDoneWork==10 or
        percDoneWork==20 or
        percDoneWork==30 or
        percDoneWork==40 or
        percDoneWork==50 or
        percDoneWork==60 or
        percDoneWork==70 or
        percDoneWork==80 or
        percDoneWork==90 or
        percDoneWork==100):
        logger.info(str(percDoneWork)+"%")
    if gui:
        sys.stdout.write("\n"+str(percDoneWork)+"%")
    else:
        sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()

def revComplement(nuc):
    return(str(Seq.Seq(nuc).reverse_complement()))

# Splits some name by _ and returns list with possible integer values
# e.g. BRAF_ex1_1_5 will give ('BRAF','ex1',int(1),int(5))
def splitNameForSorting(name):
    values=name.split('_')
    newValues=[]
    for val in values:
        try:
            newValues.append(int(val))
        except ValueError:
            newValues.append(val)
    return(newValues)

# Section of input arguments
par=argparse.ArgumentParser(description='This script constructs primers for multiplex NGS panels')
par.add_argument('--regions-file','-regions',
                 dest='regionsFile',type=str,help='file with regions for amplification in the following format:'
                 'Chromosome{Tab}Start_Position{Tab}End_Position{Tab}Amplicon_Name{Tab}\n'
                 'Desired_Multiplex_Numbers(optional){Tab}Type_Of_Primers(only left/only right/both)(optional){Tab}'
                 'Use_Whole_Region(optional)',
                 required=True)
par.add_argument('--primers-file','-primers',
                 dest='primersFile',type=str,
                 help='file with previously designed internal primers. Use this parameter, if you want only to design external primers',
                 required=False)
par.add_argument('--draft-primers','-draft',
                 dest='draftFile',type=str,
                 help='file with internal primers previously designed for part of input regions. '
                      'The program will design primers for the left regions',
                 required=False)
par.add_argument('--reference-genome','-ref',
                 dest='wholeGenomeRef',type=str,
                 help='file with INDEXED whole-genome reference sequence',
                 required=True)
par.add_argument('--adapter-for-left','-ad1',
                 dest='leftAdapter',type=str,
                 help='adapter for left primers. Use it, if you want to preserve '
                      'formation of second structures with adapter sequences (optional)',
                 required=False)
par.add_argument('--adapter-for-right','-ad2',
                 dest='rightAdapter',type=str,
                 help='adapter for right primers. Use it, if you want to preserve '
                      'formation of second structures with adapter sequences (optional)',
                 required=False)
par.add_argument('--min-amplicon-length','-minampllen',
                 dest='minAmplLen',type=int,
                 help='minimal length of amplicons. Default: 100',
                 required=False,default=100)
par.add_argument('--max-amplicon-length','-maxampllen',
                 dest='maxAmplLen',type=int,
                 help='maximal length of amplicons. Default: 110',
                 required=False,default=110)
par.add_argument('--optimal-amplicon-length','-optampllen',
                 dest='optAmplLen',type=int,
                 help='optimal length of amplicons. Default: 110',
                 required=False,default=110)
par.add_argument('--min-primer-length','-minprimerlen',
                 dest='minPrimerLen',type=int,
                 help='minimal length of primers. Default: 16',
                 required=False,default=16)
par.add_argument('--max-primer-length','-maxprimerlen',
                 dest='maxPrimerLen',type=int,
                 help='maximal length of primers. Default: 28',
                 required=False,default=28)
par.add_argument('--optimal-primer-length','-optprimerlen',
                 dest='optPrimerLen',type=int,
                 help='optimal length of primers. Default: 23',
                 required=False,default=23)
par.add_argument('--min-primer-melting-temp','-minprimermelt',
                 dest='minPrimerMelt',type=int,
                 help='minimal melting temperature of primers, degrees Celsius. Default: 60',
                 required=False,default=60)
par.add_argument('--max-primer-melting-temp','-maxprimermelt',
                 dest='maxPrimerMelt',type=int,
                 help='maximal melting temperature of primers, degrees Celsius. Default: 68',
                 required=False,default=68)
par.add_argument('--optimal-primer-melting-temp','-optprimermelt',
                 dest='optPrimerMelt',type=int,
                 help='optimal melting temperature of primers, degrees Celsius. Default: 64',
                 required=False,default=64)
par.add_argument('--min-primer-gc','-minprimergc',
                 dest='minPrimerGC',type=int,
                 help='minimal acceptable GC-content for primers. Default: 20',
                 required=False,default=20)
par.add_argument('--max-primer-gc','-maxprimergc',
                 dest='maxPrimerGC',type=int,
                 help='maximal acceptable GC-content for primers. Default: 80',
                 required=False,default=80)
par.add_argument('--optimal-primer-gc','-optprimergc',
                 dest='optPrimerGC',type=int,
                 help='optimal acceptable GC-content for primers. Default: 40',
                 required=False,default=40)
par.add_argument('--min-primer-end-gc','-minprimerendgc',
                 dest='minPrimerEndGC',type=int,
                 help="minimal acceptable number of G or C nucleotides within last 5 nucleotides of 3'-end of primers. Default: 0",
                 required=False,default=0)
par.add_argument('--max-primer-end-gc','-maxprimerendgc',
                 dest='maxPrimerEndGC',type=int,
                 help="maximal acceptable number of G or C nucleotides within last 5 nucleotides of 3'-end of primers. Default: 5",
                 required=False,default=5)
par.add_argument('--opt-primer-end-gc','-optprimerendgc',
                 dest='optPrimerEndGC',type=int,
                 help="optimal number of G or C nucleotides within last 5 nucleotides of 3'-end of primers. Default: 2",
                 required=False,default=2)
par.add_argument('--max-primer-poly-n','-maxprimerpolyn',
                 dest='maxPrimerPolyN',type=int,
                 help="maximal acceptable length of some poly-N in primers. Default: 8",
                 required=False,default=8)
par.add_argument('--max-primer-compl-end-th','-maxprimercomplendth',
                 dest='maxPrimerComplEndTh',type=int,
                 help="maximal Tm for complementarity of 3'-ends of primers. Default: 25",
                 required=False,default=25)
par.add_argument('--max-primer-compl-any-th','-maxprimercomplanyth',
                 dest='maxPrimerComplAnyTh',type=int,
                 help="maximal Tm for any complementarity of primers. Default: 35",
                 required=False,default=35)
par.add_argument('--max-primer-hairpin-th','-maxprimerhairpinth',
                 dest='maxPrimerHairpinTh',type=int,
                 help="maximal melting temperature of primer hairpin structure. Default: 40",
                 required=False,default=40)
par.add_argument('--max-primer-nonspecific','-maxprimernonspec',
                 dest='maxPrimerNonspec',type=int,
                 help="maximal number of nonspecific regions to which primer can hybridizes. Default: 10000",
                 required=False,default=10000)
par.add_argument('--max-amplicons-overlap','-maxoverlap',
                 dest='maxoverlap',type=int,
                 help='maximal length of overlap between two amplified blocks (it does not include primers). Default: 50',
                 required=False,default=50)
par.add_argument('--primers-number1','-primernum1',
                 dest='primernum1',type=int,
                 help='number of primer that user wants to get on the 1st stage. '
                      'The more this value, the more precise the choice of primers, but the longer the design time. Default: 50',
                 required=False,default=50)
par.add_argument('--auto-adjust-parameters','-autoadjust',
                 dest='autoAdjust',action='store_true',
                 help='use this parameter if you want NGS-PrimerPlex to automatically use less stringent parameters '
                      'if no primer were constructed for some region')
par.add_argument('--tries-to-get-best-combination','-tries',
                 dest='triesToGetCombination',type=int,
                 help='number of of tries to get the best primer combination. '
                      'More the value, better combination will be, '
                      'but this will take more time. Default: 10000',
                 required=False,default=10000)
par.add_argument('--return-variants-number','-returnvariantsnum',
                 dest='returnVariantsNum',type=int,
                 help='number of multiplexes variants that user wants to get after all analyses and filters. Default: 10',
                 required=False,default=10)
par.add_argument('--embedded-amplification','-embedded',
                 dest='embeddedAmpl',action='store_true',
                 help='use this parameter if you want to create NGS-panel with embedded amplification')
par.add_argument('--min-internal-primer-shift','-minprimershift',
                 dest='minPrimerShift',type=int,
                 help="minimal shift of external primer from the 3'-end of internal primer. Default: 5",
                 required=False,default=5)
par.add_argument('--opt-external-amplicon-length','-optextampllen',
                 dest='optExtAmplLen',type=int,
                 help="optimal length of the external amplicons. Default: 150",
                 required=False,default=150)
par.add_argument('--max-external-amplicon-length','-maxextampllen',
                 dest='maxExtAmplLen',type=int,
                 help="maximal length of the external amplicons. Default: 150",
                 required=False,default=150)
par.add_argument('--do-blast','-blast',
                 dest='doBlast',action='store_true',
                 help='use this parameter if you want to perform Blast-analysis of constructed primers')
par.add_argument('--substititutions-num','-subst',
                 dest='substNum',type=int,
                 help='accepted number of substitutions for searching primers in genome. Default: 2',
                 required=False,default=2)
par.add_argument('--max-nonspecific-amplicon-length','-maxnonspeclen',
                 dest='maxNonSpecLen',type=int,
                 help='maximal length of nonspecific amplicons that the program should consider. '
                      'For example, if you design primers for DNA from serum, you can set it as 150. Default: 200',
                 required=False,default=200)
par.add_argument('--snps','-snps',
                 dest='snps',action='store_true',
                 help="use this parameter if you want to check that 3'-ends of your primers do not cover any SNPs with high frequency")
par.add_argument('--dbsnp-vcf','-dbsnp',
                 dest='dbSnpVcfFile',type=str,
                 help='VCF-file (may be gzipped) with dbSNP variations',
                 required=False)
par.add_argument('--snp-freq','-freq',
                 dest='snpFreq',type=float,
                 help='minimal frequency of SNP in whole population to consider it high-frequent SNP. Default: 0.05',
                 required=False,default=0.05)
par.add_argument('--nucletide-number-to-check','-nucs',
                 dest='nucNumToCheck',type=int,
                 help='Number of nucleotides from 3`-end to check for covering SNPs. Default: None and the program will check all nucleotides',
                 required=False,default=None)
par.add_argument('--min-multiplex-dimer-dg1','-minmultdimerdg1',
                 dest='minMultDimerdG1',type=float,
                 help="minimal acceptable value of free energy of primer dimer formation "
                      "with hybridized 3'-end in one multiplex in kcal/mol. Default: -6",
                 required=False,default=-6)
par.add_argument('--min-multiplex-dimer-dg2','-minmultdimerdg2',
                 dest='minMultDimerdG2',type=float,
                 help="minimal acceptable value of free energy of primer dimer formation "
                      "in one multiplex in kcal/mol. Default: -10",
                 required=False,default=-10)
par.add_argument('--max-neighbor-intersection','-maxneighborinter',
                 dest='maxPrimerIntersection',type=float,
                 help="Maximal acceptable intersection of neighboring primers. Default: 5",
                 required=False,default=5)
par.add_argument('--threads','-th',
                 dest='threads',type=int,
                 help='number of threads. Default: 2',
                 required=False,default=2)
par.add_argument('--run-name','-run',
                 dest='runName',type=str,
                 help='name of program run. It will be used in the output file names',
                 required=False)
par.add_argument('--skip-uncovered','-skip',
                 dest='skipUndesigned',action='store_true',
                 help='use this parameter if you want to skip some targets for which primers can not be designed with defined parameters')
# Parameters for calculating thermodynamic parameters
par.add_argument('--monovalent-concentration','-mv',
                 dest='mvConc',type=int,
                 help='Concentration of monovalent cations, commonly K+ or NH4+, in mM. Default: 50',
                 required=False,default=50)
par.add_argument('--divalent-concentration','-dv',
                 dest='dvConc',type=int,
                 help='Concentration of divalent cations, commonly Mg2+, in mM. Default: 3',
                 required=False,default=3)
par.add_argument('--dntp-concentration','-dntp',
                 dest='dntpConc',type=float,
                 help='Total concentration of dNTPs. If you have each dNTP with concantration 0.2 mM, then total is 0.8 mM. Default: 0.8',
                 required=False,default=0.8)
par.add_argument('--primer-concentration','-primerconc',
                 dest='primerConc',type=int,
                 help='Concentration of each primer, in nM. Default: 250',
                 required=False,default=250)
par.add_argument('--gui','-gui',
                 dest='gui',action='store_true',
                 help='this parameter is only automatically used by GUI of the application')
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
for arg in sorted(vars(args)):
    logger.info(str(arg)+' '+str(getattr(args,arg)))
# Section of input arguments control
try:
    regionsFile=open(args.regionsFile)
except FileNotFoundError:
    print('ERROR! Input file was not found: '+args.regionsFile)
    logger.error('Input file was not found: '+args.regionsFile)
    exit(29)
if ' ' in args.regionsFile:
    print('WARNING!' +'/n'+'There is space in the directory of regions file. We recommend you to delete it to avoid errors')
if args.primersFile and args.draftFile:
    print('ERROR! You can use only primers file OR draft primers file. Leave one of this arguments')
    logger.error('ERROR! You can use only primers file OR draft primers file. Leave one of this arguments')
    exit(30)
if ((args.primersFile and
     ' ' in args.primersFile) or
    (args.draftFile and
     ' ' in args.draftFile)):
    print('WARNING!' +'/n'+'There is space in the directory of primers file. We recommend you to delete it to avoid errors')
inputDir=args.regionsFile[:args.regionsFile.rfind('/')+1]
wgref=args.wholeGenomeRef
refFa=pysam.FastaFile(args.wholeGenomeRef)
if not args.maxPrimerHairpinTh:
    args.maxPrimerHairpinTh=args.minPrimerMelt-10
if (args.embeddedAmpl or args.primersFile) and args.maxExtAmplLen<=args.maxAmplLen+2*args.minPrimerShift:
    print('ERROR! Maximal length of an extrenal amplicon should be more than maximal length of internal one plus two lengths of mininimal primer shift')
    logger.error('Maximal length of an extrenal amplicon should be more than maximal length of internal one plus two lengths of mininimal primer shift')
    exit(31)
if not os.path.exists(wgref):
    print('#'*20+'\nERROR! Whole-genome reference file does not exist:',wgref)
    logger.error('Whole-genome reference file does not exist:'+wgref)
    exit(32)
if ' ' in args.wholeGenomeRef:
    print('WARNING!' +'/n'+'There is space in the directory of reference file. We recommend you to delete it to avoid errors')
if args.snps and not args.dbSnpVcfFile:
    print('#'*20+'\nERROR! If you want to check primers for crossing SNPs, choose VCF-file with dbSNP variations (-dbsnp)!')
    logger.error('ERROR! If you want to check primers for crossing SNPs, choose VCF-file with dbSNP variations (-dbsnp)!')
    exit(33)
if args.dbSnpVcfFile and not os.path.exists(args.dbSnpVcfFile):
    print('#'*20+'\nERROR! VCF-file with dbSNP variations does not exist:',args.dbSnpVcfFile)
    logger.error('VCF-file with dbSNP variations does not exist:'+args.dbSnpVcfFile)
    exit(48)
if args.dbSnpVcfFile and ' ' in args.dbSnpVcfFile:
    print('WARNING!' +'/n'+'There is space in the directory of SNP file. We recommend you to delete it to avoid errors')

chrToChr(args,refFa)

# We make primer3 parameters for each input position
## Later we will send them into multithreading pool
print('Reading input file...')
logger.info('Reading input file...')
allRegions,regionsNames,regionsCoords,regionNameToChrom,regionNameToMultiplex,regionNameToPrimerType=readInputFile(regionsFile)

# If user also used file with primers
if args.primersFile:
    print('Reading file with primers...')
    logger.info('Reading file with primers...')
    # inputInternalPrimers and amplfiedRegions contain chromosome names
    # But regionsCoords contains chromosome numbers
    inputInternalPrimers,amplfiedRegions=readPrimersFile(args.primersFile)
    # Check that input primers cover all input regions
    chromosomesWithUncoveredRegions=[]
    for chrom,coveredRegions in amplfiedRegions.items():
        chromInt=nameToNum[chrom]
        try:
            inter=coveredRegions.intersection(regionsCoords[chromInt])
        except KeyError:
            print('ERROR! Inconsistent chromosome names in two dictionaries:')
            logger.error('Inconsistent chromosome names in two dictionaries:')
            print('regionCoords:',regionsCoords.keys())
            logger.error('regionCoords: '+str(regionsCoords.keys()))
            print('amplfiedRegions:',amplfiedRegions.keys())
            logger.error('amplfiedRegions: '+str(amplfiedRegions.keys()))
            exit(50)
        if len(inter)!=len(regionsCoords[chromInt]):
            chromosomesWithUncoveredRegions.append(chrom)
    if len(chromosomesWithUncoveredRegions)>0:
        print('ERROR! Some regions of the following chromosomes are not covered with input primers:')
        logger.error('Some regions of the following chromosomes are not covered with input primers:')
        for chrom in chromosomesWithUncoveredRegions:
            print(chrom)
            logger.error(chrom)
        exit(34)
    # If user wants to automatically sort primers pairs into multiplex reactions
    ## We create file for storing all problematic pairs of primer pairs
    if len(regionNameToMultiplex)>0:
        multiplexProblemsWB=xls.Workbook(args.primersFile[:-4]+'_amplicons_multiplex_incompatibility.xls')
        mpws=multiplexProblemsWB.add_worksheet('Internal_Primers')
        mpws.write_row(0,0,['Primers_Pair1','Primers_Pair2','Problem while joining to one multiplex'])
        colsWidth=[40,40,20]
        for k,colsWidth in enumerate(colsWidth):
            mpws.set_column(k,k,colsWidth)
    mpwsRowNum=1
    wbw=xls.Workbook(args.primersFile[:-4]+'_with_external_primers.xls')
    wsw1=wbw.add_worksheet('NGS_Primerplex_Internal_Primers')
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
        chrom=str(internalPrimers[3])
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
                        fit,problem=checkPrimersFit(internalPrimers[0:2]+[chrom,amplStart,amplEnd],
                                                    outputInternalPrimers[node],
                                                    args.minMultDimerdG1,args.minMultDimerdG2,
                                                    args.maxPrimerIntersection,[],
                                                    args.mvConc,args.dvConc,
                                                    args.dntpConc,args.primerConc,
                                                    args.leftAdapter,args.rightAdapter)
                    except KeyError:
                        print('ERROR:',outputInternalPrimers.keys())
                        print(globalMultiplexNums.nodes())
                        exit(35)
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
    print('Creating primer3 parameters for external primers design...')
    logger.info('Creating primer3 parameters for external primers design...')
    primer3Params,regionNameToChrom,amplToStartCoord=createPrimer3_parameters(allRegions,args,refFa,
                                                                              designedInternalPrimers=outputInternalPrimers)
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
            results.append(p.apply_async(runPrimer3,(regionName,inputParam,True,args)))
    outputExternalPrimers={}
    extPrimersInfo={} # This variable only for blasting designed external primers
    doneWork=0
    wholeWork=len(results)
    for res in results:
        res=res.get()
        doneWork+=1
        showPercWork(doneWork,wholeWork,args.gui)
        curRegionName,primerSeqs,primersCoords,primerTms,amplLens,amplScores,primer3Params,designOutput=res
        if primerSeqs==None:
            regionsWithoutPrimers.append(curRegionName)
            continue
        # leftPrimer,rightPrimer,amplName,chrom,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,leftPrimerTm,rightPrimerTm,leftPrimerLen,rightPrimerLen,leftGC,rightGC
        # Save all found variants of external primers. Later we will choose the best ones
        for k in range(int(len(primerSeqs)/2)):
            # We do not need to create multiplexes. We only want to surround alredy designed amplicons
            ## with external primers, selecting only the best ones
            try:
                chrom=regionNameToChrom[curRegionName]
            except:
                print('ERROR (36)!',regionNameToChrom)
                print(curRegionName)
                exit(36)
            start=amplToStartCoord[curRegionName]+primersCoords[2*k+0][0]+1
            end=amplToStartCoord[curRegionName]+primersCoords[2*k+1][0]+1
            leftGC=round(100*(primerSeqs[2*k+0].count('G')+primerSeqs[2*k+0].count('C'))/len(primerSeqs[2*k+0]),2)
            rightGC=round(100*(primerSeqs[2*k+1].count('G')+primerSeqs[2*k+1].count('C'))/len(primerSeqs[2*k+1]),2)
            chromName=chrom
            chromInt=nameToNum[chrom]
            if start-1-100<1:
                extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                 chromName,1,end+100)
            elif end+100>refFa.lengths[refFa.references.index(chromName)]:
                extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                 chromName,start-1-100,refFa.lengths[refFa.references.index(chromName)])
            else:
                extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                 chromName,start-1-100,end+100)
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
            print('\n # WARNING! For',list(regionsWithoutPrimersCounter.values()).count(3),'amplicon(s) external primers were not designed! Try less stringent parameters.')
            logger.warn(' # WARNING! For '+str(list(regionsWithoutPrimersCounter.values()).count(3))+' amplicon(s) external primers were not designed! Try less stringent parameters.')
            for regionsWithoutPrimer,value in sorted(regionsWithoutPrimersCounter.items(),key=itemgetter(1),reverse=True):
                if value<3: break
                print('   '+regionsWithoutPrimer)
                logger.info('   '+regionsWithoutPrimer)
    print('\n # Number of designed external primers: '+str(externalPrimersNum))
    logger.info(' # Number of designed external primers: '+str(externalPrimersNum))
    p.close()
    p.join()
    writeDraftPrimers(extPrimersInfo,
                      args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers.xls',
                      external=True)
    # Analyzing external primers specificity        
    if args.doBlast:
        print('\nAnalyzing external primers for their specificity...')
        logger.info('Analyzing external primers for their specificity...')
        specificPrimers,primersNonSpecRegionsByChrs=checkPrimersSpecificity(args.regionsFile[:-4],extPrimersInfo,args.wholeGenomeRef,
                                                                            args.runName,refFa,args.substNum,args.threads,args.gui,
                                                                            args.maxNonSpecLen,args.maxPrimerNonspec,True,str(i+1))
        print(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
        logger.info(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
        # Write primers that left after filtering by specificity
        writeDraftPrimers(extPrimersInfo,
                          args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_specificity.xls',
                          goodPrimers=specificPrimers,
                          external=True)
    # Check external primers for covering high-frequent SNPs
    if args.snps:
        print('Analyzing external primers for covering high-frequent SNPs...')
        logger.info('Analyzing external primers for covering high-frequent SNPs...')
        primerPairsNonCoveringSNPs,primersCoveringSNPs=analyzePrimersForCrossingSNP(extPrimersInfo,
                                                                                    args.threads,
                                                                                    args.dbSnpVcfFile,
                                                                                    args.nucNumToCheck,
                                                                                    args.gui)
        print("\n # Number of primers covering high-frequent SNPs: "+str(len(primersCoveringSNPs))+'. They will be removed.')
        logger.info(" # Number of primers covering high-frequent SNPs: "+str(len(primersCoveringSNPs))+'. They will be removed.')
        # Write primers that left after filtering by SNPs
        writeDraftPrimers(extPrimersInfo,
                          args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_SNPs.xls',
                          goodPrimers=primerPairsNonCoveringSNPs,
                          external=True)
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
        unspecificPrimers=getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,args.maxNonSpecLen,
                                                                  args.threads,args.gui)
        print(' # Number of external primer pairs that form unspecific product:',len(unspecificPrimers))
        logger.info(' # Number of external primer pairs that form unspecific product: '+str(len(unspecificPrimers)))
    elif not args.doBlast:
        unspecificPrimers=[]
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
                fit,problem=checkPrimersFit(extPrimers[0:2]+extPrimers[3:6],
                                            outputExternalPrimers[node],
                                            args.minMultDimerdG1,args.minMultDimerdG2,
                                            args.maxPrimerIntersection,unspecificPrimers,
                                            args.mvConc,args.dvConc,
                                            args.dntpConc,args.primerConc)
                if not fit:
                    mpws.write_row(mpwsRowNum,0,[','.join(extPrimers[0:2]),','.join(outputExternalPrimers[node][0:2]),problem])
                    mpwsRowNum+=1
                if not fit and globalMultiplexNums.has_edge(node,regionName):
                    globalMultiplexNums.remove_edge(node,regionName)
        multiplexProblemsWB.close()
    for chrom,coords in sorted(primersForOutput.items(),
                               key=lambda item:(nameToNum[item[0]],
                                                item[1])):
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
                    print('ERROR #38 ',amplNameToRowNum[regionName],primer[:-1])
                    logger.error('#38 '+str(amplNameToRowNum[regionName])+'\n'+str(primer[:-1]))
                    exit(38)
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
            cls=makeFinalMultiplexes(localMultiplexNums,[],len(mults),True)
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
                leftNodesGraph=deepcopy(localMultiplexNums)
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
                exit(39)
        for k,multiplex in enumerate(multiplexes):
            for ampl in multiplex:
                wsw2.write(amplNameToRowNum[ampl],18,k+1)
                wsw1.write(amplNameToRowNum[ampl],17,k+1)
    print()
    logger.info('\n')
    wbw.close()
else:
    # If user use as input draft primers
    if args.draftFile:
        print('Reading file with draft primers...')
        logger.info('Reading file with draft primers...')
        # Read file with draft primers
        primersInfo,primersInfoByChrom=readDraftPrimers(args.draftFile)
        print('Number of primer pairs from draft file: '+str(len(primersInfo)))
        logger.info('Number of primer pairs from draft file: '+str(len(primersInfo)))
##        if not args.skipUndesigned:
        # Get regions that uncovered by already designed primers
        print('Getting positions that uncovered by draft primers...')
        logger.info('Getting positions that uncovered by draft primers...')
        uncoveredRegions,amplNames,primersToAmplNames=getRegionsUncoveredByDraftPrimers(allRegions,primersInfoByChrom)
        if len(uncoveredRegions)>0:
            # Go through all regions sorted by chromosome and coordinate of start
            print('Creating input parameters for primer3...')
            logger.info('Creating input parameters for primer3...')
            primer3Params=createPrimer3_parameters(uncoveredRegions,args,refFa,
                                                   regionNameToPrimerType=regionNameToPrimerType)
            # Construct primers for each created set of parameters
            print('Constructing primers...')
            logger.info('Constructing primers...')
            # primer3Params,regionNameToChrom,args,regionsCoords=None,allRegions=None,primersInfo=None,primersInfoByChrom=None,amplNames=None,primersToAmplNames=None
            primersInfo,primersInfoByChrom,amplNames,primersToAmplNames,regionsCoords,regionNameToChrom,uncoveredRegions=constructInternalPrimers(primer3Params,regionNameToChrom,args,regionsCoords,allRegions,
                                                                                                                                                  primersInfo,primersInfoByChrom,amplNames,primersToAmplNames)
            print('Total number of primer pairs: '+str(len(primersInfo)))
            logger.info('Total number of primer pairs: '+str(len(primersInfo)))
    else:                            
        # Go through all regions sorted by chromosome and coordinate of start
        print('Creating input parameters for primer3...')
        logger.info('Creating input parameters for primer3...')
        primer3Params=createPrimer3_parameters(allRegions,args,refFa,regionNameToPrimerType=regionNameToPrimerType)

        # Construct primers for each created set of parameters
        print('Constructing primers...')
        logger.info('Constructing primers...')
        primersInfo,primersInfoByChrom,amplNames,primersToAmplNames,regionsCoords,regionNameToChrom,allRegions=constructInternalPrimers(primer3Params,regionNameToChrom,args,regionsCoords,allRegions)

    # Check all primers with BWA
    if args.doBlast:
        print('Analyzing primers for their specificity...')
        logger.info('Analyzing primers for their specificity...')
        specificPrimers,primersNonSpecRegionsByChrs=checkPrimersSpecificity(args.regionsFile[:-4],primersInfo,args.wholeGenomeRef,
                                                                            args.runName,refFa,args.substNum,args.threads,args.gui,
                                                                            args.maxNonSpecLen,args.maxPrimerNonspec,False,'')    
        print(' # Number of specific primer pairs:',len(specificPrimers))
        logger.info(' # Number of specific primer pairs: '+str(len(specificPrimers)))
        # Now we need to remove all unspecific primers from constructed primer pairs
        print(' Removing unspecific primer pairs...')
        logger.info(' Removing unspecific primer pairs...')
        primersInfoByChrom,primersInfo,amplNames=removeBadPrimerPairs(primersInfoByChrom,primersInfo,specificPrimers,primersToAmplNames,amplNames)
        # Write primers that left after filtering by specificity
        writeDraftPrimers(primersInfo,args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_specificity.xls',specificPrimers)
        # If we filtered primers out by specificity, we need to check that all input regions are still covered
        amplNames,allRegions,regionNameToChrom,regionsCoords=checkThatAllInputRegionsCovered(amplNames,allRegions,regionNameToChrom,regionsCoords,'by specificity',args.skipUndesigned)
        # If user wants to automatically sort amplicons by multiplexes
        if len(regionNameToMultiplex)>0:
            # We save primers from different pairs that form unspecific amplicons (for sorting primer pairs by multiplexs later)
            print(' Searching for nonspecific amplicons that are formed by primers from different primer pairs...')
            logger.info(' Searching for nonspecific amplicons that are formed by primers from different primer pairs...')
            unspecificPrimers=getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,args.maxNonSpecLen,
                                                                      args.threads,args.gui)
            print('\n # Number of primer pairs that form unspecific product:',len(unspecificPrimers))
            logger.info(' # Number of primer pairs that form unspecific product: '+str(len(unspecificPrimers)))
##            writeUnspecificPrimers(primersInfo,
##                                   args.regionsFile[:-4]+'_NGS_primerplex_unspecific_products.xls',
##                                   unspecificPrimers)
        else:
            unspecificPrimers=[]
    else:
        unspecificPrimers=[]
    # Check primers for covering high-frequent SNPs
    if args.snps:
        print('Analyzing primers for covering high-frequent SNPs...')
        logger.info('Analyzing primers for covering high-frequent SNPs...')
        primerPairsNonCoveringSNPs,primersCoveringSNPs=analyzePrimersForCrossingSNP(primersInfo,
                                                                                    args.threads,
                                                                                    args.dbSnpVcfFile,
                                                                                    args.nucNumToCheck,
                                                                                    args.gui)
        if len(primersCoveringSNPs)>0:
            print(" Removing primer pairs covering high-frequent SNPs...")
            logger.info(" Removing primer pairs covering high-frequent SNPs...")
            primersInfoByChrom,primersInfo,amplNames=removeBadPrimerPairs(primersInfoByChrom,primersInfo,primerPairsNonCoveringSNPs,primersToAmplNames,amplNames)
            # Write primers that left after filtering by covering SNPs
            writeDraftPrimers(primersInfo,args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_SNPs.xls',primerPairsNonCoveringSNPs)
            # If we filtered primers out that cover SNPs, we need to check that all input regions are still covered
            amplNames,allRegions,regionNameToChrom,regionsCoords=checkThatAllInputRegionsCovered(amplNames,allRegions,regionNameToChrom,regionsCoords,'that cover SNPs',args.skipUndesigned)
        else:
            writeDraftPrimers(primersInfo,args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_SNPs.xls',primerPairsNonCoveringSNPs)

    # Now we need to group primers into continuous amplified blocks
    # amplified block contains all multiplexes for ditinct joined regions
    ## [[joined_regions_1],[joined_region_2],[joined_region_3]...]
    ### Two different joined regions may be located on one or two chromosomes
    allRegionsAmplifiedBlocks={}
    pPolyN=re.compile('(A+|C+|T+|G+)')
    # We go through all chromosomes
    print('Joining primer pairs to amplified blocks...')
    logger.info('Joining primer pairs to amplified blocks...')
    p=ThreadPool(args.threads)
    results=[]
    for chromInt in sorted(regionsCoords.keys()):
        chrom=numToName[chromInt]
        results.append(p.apply_async(joinAmpliconsToBlocks,(sorted(regionsCoords[chromInt]),
                                                            primersInfoByChrom[chromInt],
                                                            args.maxAmplLen,chromInt,
                                                            args.returnVariantsNum)))
    doneWork=0
    wholeWork=len(results)
    showPercWork(doneWork,wholeWork,args.gui)
    for res in results:
        chromInt,finalShortestPaths=res.get()
        allRegionsAmplifiedBlocks[chromInt]=finalShortestPaths
        doneWork+=1
        showPercWork(doneWork,wholeWork,args.gui)
    # We need to show user statistics of amplified blocks
    totalAmplicons=0
    totalMultiplexVariants=1
    blockNum=0
    blockToChromNum={}
    # One chromosome may contain severel "amplified blocks"
    # Later we can combine primer pairs from different amplified blocks in different variants
    for chromInt,blocks in sorted(allRegionsAmplifiedBlocks.items()):
        chrom=numToName[chromInt]
        print('\n # Total number of amplified blocks on chromosome',str(chrom)+':',len(allRegionsAmplifiedBlocks[chromInt]))
        logger.info(' # Total number of amplified blocks on chromosome '+str(chrom)+': '+str(len(allRegionsAmplifiedBlocks[chromInt])))
        for i,block in enumerate(blocks):
            blockToChromNum[blockNum]=[chromInt,i]
            blockNum+=1
            totalMultiplexVariants*=len(block)
            for mult in block:
                totalAmplicons+=len(mult)
    print(' # Total number of amplified blocks:',blockNum)
    print(' # Total number of amplicons:',totalAmplicons)
    logger.info(' # Total number of amplified blocks: '+str(blockNum))
    logger.info(' # Total number of amplicons: '+str(totalAmplicons))
    if totalMultiplexVariants>10**6:
        print(' # Total number of multiplex variants: 10^'+str(int(round(math.log10(totalMultiplexVariants),0))))
        logger.info(' # Total number of multiplex variants: 10^'+str(int(round(math.log10(totalMultiplexVariants),0))))
    else:
        print(' # Total number of multiplex variants: '+str(totalMultiplexVariants))
        logger.info(' # Total number of multiplex variants: '+str(totalMultiplexVariants))
    # Get best combinations of amplified block variants
    combinations=getBestPrimerCombinations(allRegionsAmplifiedBlocks,
                                           primersInfo,
                                           unspecificPrimers,
                                           args.returnVariantsNum,
                                           args.triesToGetCombination,
                                           totalMultiplexVariants,
                                           args.gui)
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
        # Chromosome in combination is chromosome number but in string format
        for chromIntStr,coords in sorted(comb.items(),
                                         key=lambda item:int(item[0])):
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
                chromName=str(numToName[int(chromIntStr)])
                if amplStart-1-100<1:
                    extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                     chromName,1,amplEnd+100)
                elif amplEnd+100>refFa.lengths[refFa.references.index(chromName)]:
                    extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                     chromName,amplStart-1-100,refFa.lengths[refFa.references.index(chromName)])
                else:
                    extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                     chromName,amplStart-1-100,amplEnd+100)
##                out=pysam.faidx(wgref,chromName+':'+str(amplStart-1-100)+'-'+str(amplEnd+100))
##                lines=out.split('\n')
##                if lines[1]=='':
##                    print('ERROR! Extracted sequence has no length')
##                    print(chromName,amplStart-1-100,amplEnd+100)
##                    exit(40)
##                extendedAmplSeq=''.join(lines[1:-1]).upper()
                rInternalFile.write('\n'.join(['>'+primersToAmplNames['_'.join(c)][0],extendedAmplSeq])+'\n')
                amplName=primersToAmplNames['_'.join(c)][0]
                regionName=amplName[:amplName.rfind('_')]
                if len(regionNameToMultiplex)>0:
                    wsw1.write_row(rowNum,0,[rowNum,c[0],c[1],amplName,chromName,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,
                                            leftPrimerTm,rightPrimerTm,len(c[0]),len(c[1]),leftGC,rightGC,','.join(regionNameToMultiplex[regionName])])
                else:
                    wsw1.write_row(rowNum,0,[rowNum,c[0],c[1],amplName,chromName,amplStart,amplEnd,amplLen,amplBlockStart,amplBlockEnd,
                                            leftPrimerTm,rightPrimerTm,len(c[0]),len(c[1]),leftGC,rightGC])
                outputInternalPrimers[primersToAmplNames['_'.join(c)][0]]=[c[0],c[1],amplName,chromName,
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
                        currentGlobalMultiplexNums=deepcopy(globalMultiplexNums)
                        for node in currentGlobalMultiplexNums.nodes():
                            try:
                                fit,problem=checkPrimersFit(c+[chromName,amplStart,amplEnd],
                                                            outputInternalPrimers[node],
                                                            args.minMultDimerdG1,args.minMultDimerdG2,
                                                            args.maxPrimerIntersection,unspecificPrimers,
                                                            args.mvConc,args.dvConc,
                                                            args.dntpConc,args.primerConc,
                                                            args.leftAdapter,args.rightAdapter)
                            except KeyError:
                                print('ERROR:',outputInternalPrimers.keys())
                                print(globalMultiplexNums.nodes())
                                exit(41)
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
            primer3Params={}
            if args.draftFile:
                print('Reading file with draft primers...')
                logger.info('Reading file with draft primers...')
                # Read file with draft primers
                extPrimersInfo,extPrimersInfoByChrom=readDraftPrimers(args.draftFile,external=True)
                print('Number of external primer pairs from draft file: '+str(len(extPrimersInfo)))
                logger.info('Number of external primer pairs from draft file: '+str(len(extPrimersInfo)))
                # Get regions that uncovered by already designed primers
                print('Getting positions that uncovered by draft external primers...')
                logger.info('Getting positions that uncovered by draft external primers...')
                # primersInfo,primersInfoByChrom,outputInternalPrimers,args
                output=getRegionsUncoveredByDraftExternalPrimers(extPrimersInfo,extPrimersInfoByChrom,
                                                                 outputInternalPrimers,args,refFa)
                outputExternalPrimers,uncoveredInternalPrimers,regionNameToChrom,amplToStartCoord=output
                print('Number of internal amplicons that do not have external primers: '+str(len(uncoveredInternalPrimers)))
                logger.info('Number of internal amplicons that do not have external primers: '+str(len(uncoveredInternalPrimers)))
                if len(uncoveredInternalPrimers)>0:
                    # Go through all regions sorted by chromosome and coordinate of start
                    print('Creating input parameters for primer3...')
                    logger.info('Creating input parameters for primer3...')
                    primer3Params,regionNameToChrom,amplToStartCoord=createPrimer3_parameters(uncoveredRegions,args,refFa,
                                                                                        designedInternalPrimers=uncoveredInternalPrimers,
                                                                                        regionNameToPrimerType=regionNameToPrimerType)
            else:
                print('\nCreating primer3 parameters for external primers design, combination variant '+str(i+1)+'...')
                logger.info('Creating primer3 parameters for external primers design, combination variant '+str(i+1)+'...')
                primer3Params,regionNameToChrom,amplToStartCoord=createPrimer3_parameters(allRegions,args,refFa,
                                                                                    designedInternalPrimers=outputInternalPrimers,
                                                                                    regionNameToPrimerType=regionNameToPrimerType,
                                                                                    regionNameToChrom=regionNameToChrom)
                outputExternalPrimers={}
                extPrimersInfo={} # This variable only for blasting designed external primers and for draft-files

##            createExternalPrimers(primer3Params,regionNameToChrom,amplToStartCoord,wgref,args,wbw,colsWidth2)
            externalPrimersNum=0
            regionsWithoutPrimers=[]
            primersForOutput={}
            wsw2=wbw.add_worksheet('NGS_Primerplex_External_Primers')
            wsw2.write_row(0,0,['#','Left_Primer_Seq','Right_Primer_Seq','Amplicon_Name','Chrom','Amplicon_Start','Amplicon_End','Amplicon_Length',
                               'Amplified_Block_Start','Amplified_Block_End','Left_Primer_Tm','Right_Primer_Tm','Left_Primer_Length','Right_Primer_Length','Left_GC','Right_GC',"Left_3'-shift","Right_3'-shift",'Multiplex','Specificity'])
            for k,colsWidth in enumerate(colsWidth2):
                wsw2.set_column(k,k,colsWidth)
            if len(primer3Params)>0:
                p=ThreadPool(args.threads)                
                print('Constructing external primers...')
                logger.info('Constructing external primers...')
                results=[]
                for regionName,inputParams in primer3Params.items():
                    for inputParam in inputParams:
                        results.append(p.apply_async(runPrimer3,(regionName,inputParam,True,args)))
                doneWork=0
                wholeWork=len(results)
                primerDesignExplains={}
                for res in results:
                    doneWork+=1
                    showPercWork(doneWork,wholeWork,args.gui)
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
                            chrom=regionNameToChrom[curRegionName]
                        except:
                            print('ERROR!',regionNameToChrom)
                            print(curRegionName)
                            exit(42)
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
                        chromName=chrom
                        if start-1-100<1:
                            extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                             chromName,1,end+100)
                        elif end+100>refFa.lengths[refFa.references.index(chromName)]:
                            extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                             chromName,start-1-100,refFa.lengths[refFa.references.index(chromName)])
                        else:
                            extendedAmplSeq=extractGenomeSeq(refFa,args.wholeGenomeRef,
                                                             chromName,start-1-100,end+100)
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
                        print(' # WARNING! For',list(regionsWithoutPrimersCounter.values()).count(3),'amplicon(s) external primers were not designed! Try less stringent parameters.')
                        logger.warn(' # WARNING! For '+str(list(regionsWithoutPrimersCounter.values()).count(3))+' amplicon(s) external primers were not designed! Try less stringent parameters.')
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
            writeDraftPrimers(extPrimersInfo,
                              args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers.xls',
                              external=True)
            # Analyzing external primers specificity        
            if args.doBlast:
                print('\nAnalyzing external primers for their specificity...')
                logger.info('Analyzing external primers for their specificity...')
                specificPrimers,primersNonSpecRegionsByChrs=checkPrimersSpecificity(args.regionsFile[:-4],extPrimersInfo,args.wholeGenomeRef,
                                                                                    args.runName,refFa,args.substNum,args.threads,args.gui,
                                                                                    args.maxNonSpecLen,args.maxPrimerNonspec,True,str(i+1))
                print(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
                logger.info(' # Number of specific external primer pairs: '+str(len(specificPrimers))+'. Unspecific pairs will be removed.')
                # Write primers that left after filtering by specificity
                writeDraftPrimers(extPrimersInfo,
                                  args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_specificity.xls',
                                  goodPrimers=specificPrimers,
                                  external=True)
            # Check external primers for covering high-frequent SNPs
            if args.snps:
                print('Analyzing external primers for covering high-frequent SNPs...')
                logger.info('Analyzing external primers for covering high-frequent SNPs...')
                primerPairsNonCoveringSNPs,primersCoveringSNPs=analyzePrimersForCrossingSNP(extPrimersInfo,
                                                                                            args.threads,
                                                                                            args.dbSnpVcfFile,
                                                                                            args.nucNumToCheck,
                                                                                            args.gui)
                print("\n # Number of primers covering high-frequent SNPs: "+str(len(primersCoveringSNPs))+'. They will be removed.')
                logger.info(" # Number of primers covering high-frequent SNPs: "+str(len(primersCoveringSNPs))+'. They will be removed.')
                # Write primers that left after filtering by SNPs
                writeDraftPrimers(extPrimersInfo,
                                  args.regionsFile[:-4]+'_NGS_primerplex_all_draft_primers_after_SNPs.xls',
                                  goodPrimers=primerPairsNonCoveringSNPs,
                                  external=True)
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
                unspecificPrimers=getPrimerPairsThatFormUnspecificProduct(primersNonSpecRegionsByChrs,args.maxNonSpecLen,
                                                                          args.threads,args.gui)
                print(' # Number of external primer pairs that form unspecific product:',len(unspecificPrimers))
                logger.info(' # Number of external primer pairs that form unspecific product: '+str(len(unspecificPrimers)))
            elif not args.doBlast:
                unspecificPrimers=[]
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
                        fit,problem=checkPrimersFit(extPrimers[0:2]+extPrimers[3:6],
                                                    outputExternalPrimers[node],
                                                    args.minMultDimerdG1,args.minMultDimerdG2,
                                                    args.maxPrimerIntersection,unspecificPrimers,
                                                    args.mvConc,args.dvConc,
                                                    args.dntpConc,args.primerConc)
                        if not fit:
                            mpws.write_row(mpwsRowNum,0,[','.join(extPrimers[0:2]),','.join(outputExternalPrimers[node][0:2]),problem])
                            mpwsRowNum+=1
                        if not fit and globalMultiplexNums.has_edge(node,regionName):
                            globalMultiplexNums.remove_edge(node,regionName)
                multiplexProblemsWB.close()
            rowNum=1
            rExternalFile=open(inputFileBase+'_NGS_primerplex'+runName+'_primers_combination_'+str(i+1)+'_external_amplicons.fa','w')
            for chrom,coords in sorted(primersForOutput.items(),
                                       key=lambda item:(nameToNum[item[0]],
                                                        item[1])):
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
                            exit(44)
            rExternalFile.close()
            if len(regionNameToMultiplex)>0:
                # Contains amplicons for each of multiplex
                multiplexes=[]
                # Contains amplicons that were not sorted to any multiplex
                leftUnsortedAmpls=[]
                # globalMultiplexesContainer different groups of multiplex numbers
                # e.g. it can be 1_2_3 for some genes and 4_5_6 for other genes
                for k,(key,containerNodes) in enumerate(sorted(globalMultiplexesContainer.items())):
                    print('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
                    logger.info('Sorting amplicons to multiplexes '+str(key)+' (contain '+str(len(containerNodes+leftUnsortedAmpls))+' amplicons)...')
                    mults=key.split('_')
                    # Create graph with all left unsorted amplicons
                    localMultiplexNums=nx.Graph()
                    localMultiplexNums=globalMultiplexNums.subgraph(containerNodes+leftUnsortedAmpls)
                    # Contains amplicons that were sorted to current multiplex
                    cls=[]
                    cls=makeFinalMultiplexes(localMultiplexNums,[],len(mults),True)
                    multiplexes.extend(cls)
                    # Calculate how many amplicons were sorted to current multiplex
                    sumLen=sum([len(x) for x in cls])
                    allSortedAmpls=[]
                    for x in cls:
                        allSortedAmpls.extend(x)
                    # If all amplicons were sorted to multiplexes
                    if sumLen==len(localMultiplexNums):
                        print(' # Number of multiplexes:',len(cls))
                        logger.info(' # Number of multiplexes: '+str(len(cls)))
                        for z,cl in enumerate(cls):
                            print('  # Number of amplicons in multiplex',mults[z],len(cl))
                            logger.info('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
                    # If not all amplicons were sorted to multiplexes
                    elif sum([len(x) for x in cls])<len(localMultiplexNums):
                        print('Number of designed multiplexes:',len(cls))
                        logger.warn(' # Number of designed multiplexes: '+str(len(cls)))
                        for z,cl in enumerate(cls):
                            print('Number of amplicons in multiplex',mults[z],len(cl))
                            logger.warn('  # Number of amplicons in multiplex '+mults[z]+': '+str(len(cl)))
                        leftNodesGraph=deepcopy(localMultiplexNums)
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
                        exit(45)
                for k,multiplex in enumerate(multiplexes):
                    for ampl in multiplex:
                        wsw2.write(amplNameToRowNum[ampl],18,k+1)
                        wsw1.write(amplNameToRowNum[ampl],17,k+1)
        wbw.close()
        print()
        logger.info('\n')
print('NGS-PrimerPlex finished!')
logger.info('NGS-PrimerPlex finished!')

# TODO:
# make function that writes list of unspecific products into file
# make function that reads list of unspecific products from file
