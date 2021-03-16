# This script checks primers efficiency with all possible approaches

import sys
import argparse
import logging
import xlrd
import primer3
import pysam
import xlsxwriter as xls
from multiprocessing.pool import ThreadPool
from xlrd.biffh import XLRDError
from Bio.SeqUtils import GC
from Bio.Seq import Seq
# Own modules
from primer_non_target_amplicons import getNonTargetAmplicons
from primer_non_targets import checkPrimersNonTargets
from work_percentage_process import showPercWork

# Functions
def revComplement(nuc):
    return(str(Seq(nuc).reverse_complement()))

def checkHairpinEnd(primer):
    primerEnd3_rc=revComplement(primer[-4:])
    primerHairpin=round(primer3.calcHairpin(primer,
                                            mv_conc=args.mvConc,
                                            dv_conc=args.dvConc).dg/1000,2)
    if primerEnd3_rc in primer[1:-3-4] and primerHairpin<-1:
        return(str(primerHairpin)+'*')
    return(primerHairpin)

def checkDimerEnd(primer1,primer2=None):
    if not primer2:
        primer2=primer1
    primersDimer=round(primer3.calcHeterodimer(primer1,
                                               primer2,
                                               mv_conc=args.mvConc,
                                               dv_conc=args.dvConc).dg/1000,2)
    primer1End_rc=revComplement(primer1[-5:])
    primer2End_rc=revComplement(primer2[-5:])
    if (primersDimer<-3 and
        (primer1End_rc in primer2[1:] or
         primer2End_rc in primer1[1:])):
        return(str(primersDimer)+'*')
    return(primersDimer)

def extractGenomeSeq(refFa,chrom,start,end):
    attempts=0
    seq=None
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

def getAmpliconGC(wholeGenomeRef,chrom,start,end):
    refFa=pysam.FastaFile(wholeGenomeRef)
    seq=extractGenomeSeq(refFa,chrom,start,end)
    gcContent=round((seq.count('C')+seq.count('G'))*100/len(seq),2)
    return(gcContent)

def getAmpliconAG(wholeGenomeRef,chrom,start,end):
    refFa=pysam.FastaFile(wholeGenomeRef)
    seq=extractGenomeSeq(refFa,chrom,start,end)
    agContent=round((seq.count('A')+seq.count('G'))*100/len(seq),2)
    return(agContent)

def getAmpliconPalindromes(wholeGenomeRef,chrom,start,end,kmerSize):
    refFa=pysam.FastaFile(wholeGenomeRef)
    seq=extractGenomeSeq(refFa,chrom,start,end)
    kmers=[]
    for i in range(len(seq)-kmerSize+1):
        kmers.append(seq[i:i+kmerSize])
    palindromesNum=0
    for i,kmer in enumerate(kmers):
        kmer_rc=revComplement(kmer)
        if kmer_rc in kmers[i+1:]:
            palindromesNum+=kmers[i+1:].count(kmer_rc)
    return(palindromesNum)

# Reading arguments
par=argparse.ArgumentParser(description='This tool checks primers efficiency with all possible approaches')
par.add_argument('--input-file','-in',
                 dest='inputFile',type=str,
                 help='Input tabulated- or excel-file '
                      'that contains sequences of primers. '
                      'One line - one amplicon primers: '
                      'Amplicon_Number, Primer1, Primer2, Amplicon_Name,'
                      'Chromosome, Amplicon_Start, Amplicon_End.'
                      'The last three columns are not necessary',
                 required=True)
par.add_argument('--whole-genome-reference','-wgref',
                 dest='wholeGenomeRef',type=str,
                 help='file with INDEXED whole-genome reference sequence.'
                      'It is necessary only if you want to blast your primers',
                 required=False)
par.add_argument('--max-non-target-amplicon-length','-maxnontargetampllen',
                 dest='maxNonTargetLen',type=int,
                 help='maximal length of non-target amplicons '
                      'that the program should consider. '
                      'For example, if you design primers for DNA from serum,'
                      'you can set it as 150. Default: 200',
                 required=False,default=200)
par.add_argument('--max-primer-nontargets','-maxprimernontargets',
                 dest='maxPrimerNonTargets',type=int,
                 help='maximal number of nontarget regions '
                      'to which primer can can be mapped '
                      'and the program will still search for '
                      'non-target amplicons formed. Default: 1000',
                 required=False,default=1000)
par.add_argument('--substititutions-num','-subst',
                 dest='substNum',type=int,
                 help='accepted number of substitutions for searching'
                      'primers in genome. Default: 2',
                 required=False,default=2)
par.add_argument('--blast-num','-bnum',
                 dest='blastNum',type=int,
                 help='number of BWA runs to find all primer targets. It is recommended to increase for genes with many pseodugenes (e.g. CHEK2 or PMS2). Default: 1',
                 required=False,default=1)
par.add_argument('--threads','-th',dest='threads',type=int,
                 help='number of threads. Default: 2',
                 required=False,default=2)
par.add_argument('--output-type','-type',dest='outType',type=str,
                 help='Type of output (xls or tsv). Default=xls',
                 default='xls')
# Parameters for calculating thermodynamic parameters
par.add_argument('--monovalent-concentration','-mv',
                 dest='mvConc',type=int,
                 help='Concentration of monovalent cations,'
                      'commonly K+ or NH4+, in mM. Default: 50',
                 required=False,default=50)
par.add_argument('--divalent-concentration','-dv',
                 dest='dvConc',type=int,
                 help='Concentration of divalent cations,'
                      'commonly Mg2+, in mM. Default: 3',
                 required=False,default=3)
par.add_argument('--dntp-concentration','-dntp',
                 dest='dntpConc',type=float,
                 help='Total concentration of dNTPs. '
                      'If you have each dNTP with concantration 0.2 mM,'
                      'then total is 0.8 mM. Default: 0.8',
                 required=False,default=0.8)
par.add_argument('--primer-concentration','-primerconc',
                 dest='primerConc',type=int,
                 help='Concentration of each primer, in nM. Default: 300',
                 required=False,default=300)
global args
args=par.parse_args()

# Set base for output files
inputFileBase=args.inputFile[:-4]

# Set logging
logger=logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler=logging.FileHandler(inputFileBase+'_checkPrimers.log')
handler.setLevel(logging.INFO)
formatter=logging.Formatter(' - '.join(['%(asctime)s',
                                        '%(name)s',
                                        '%(levelname)s'
                                        '%(message)s']))
handler.setFormatter(formatter)
logger.addHandler(handler)
# Write command to the log
logger.info('The command was:\n'+' '.join(sys.argv))
if args.outType not in ['xls','tsv']:
    print('ERROR! Unknown type of output! It should be tsv or xls')
    exit(1)

# Reading input
primerSets={}
primerSetOrder=[]
primerInfo={}
try:
    wb=xlrd.open_workbook(args.inputFile)
    ws=wb.sheet_by_index(0)
    for i in range(ws.nrows):
        row=ws.row_values(i)
        if 'Name' in row or i==0:
            continue
        if row[0]=='':
            break
        if len(row)>4 and 'chr' not in str(row[4]):
            try:
                row[4]='chr'+str(int(row[4]))
            except ValueError:
                pass
        if (int(row[0]) in primerSets.keys() or
            int(row[0]) in primerInfo.keys()):
            print('ERROR! Numbers of amplicons are repeated:',int(row[0]))
            exit(1)
        if len(row)<4:
            row.append(row[0])
        elif row[3]=='':
            row[3]=row[0]
        primerSets[int(row[0])]='_'.join(row[1:3])
        primerInfo[int(row[0])]=[str(int(row[0])),*row[1:]]
        if int(row[0]) not in primerSetOrder:
            primerSetOrder.append(int(row[0]))
        else:
            print('ERROR! Names of amplicons are repeated:',int(row[0]))
            exit(1)
except XLRDError:
    file=open(args.inputFile)
    for string in file:
        if 'Name' in string: continue
        row=string.replace('\n','').split('\t')
        if 'chr' not in str(row[4]):
            row[4]='chr'+str(int(row[4]))
        if (int(row[0]) in primerSets.keys() or
            int(row[0]) in primerInfo.keys()):
            print('ERROR! Numbers of amplicons are repeated:',int(row[0]))
            exit(1)
        primerSets[int(row[0])]='_'.join(row[1:3])
        primerInfo[int(row[0])]=row
        if int(row[0]) not in primerSetOrder:
            primerSetOrder.append(int(row[0]))
        else:
            print('ERROR (1)! Names of amplicons are repeated:',int(row[0]))
            exit(1)
    file.close()

# Checking primers and writing result to output
if args.outType=='xls':
    wbw=xls.Workbook(inputFileBase+'_checkPrimers_result.xls')
    wsw=wbw.add_worksheet()
elif args.outType=='tsv':
    rFile=open(inputFileBase+'_checkPrimers_result.csv','w')
colNames=['Amplicon_Number','Primer1','Primer2',
          'Amplicon_Name','Chromosome','Start','End','Amplicon_Length',
          'Amplified_Block_Start','Amplified_Block_End',
          'Length1','Length2','Tm1','Tm2','GC1','GC2',
          'end3_GC1','end3_GC2','Amplicon_GC','Amplicon_AG',
          'Amplicon_2mer_Palindromes',
          'Amplicon_3mer_Palindromes','Amplicon_4mer_Palindromes',
          'Amplicon_5mer_Palindromes','Amplicon_6mer_Palindromes',
          'Amplicon_7mer_Palindromes','Amplicon_8mer_Palindromes',
          'Hairpin1','Hairpin2','Homodimer1','Homodimer2','Heterodimer']
if args.wholeGenomeRef:
    colNames+=['Nonspecific_Regions1',
               'Nonspecific_Regions2','Nonspecific_Amplicons']
colWidths=[7,30,30,10,10,10,10,10,
           10,10,
           5,5,5,5,5,5,
           9,9,11,11,
           20,
           20,20,
           20,20,
           20,20,
           8,8,13,13,13,
           25,18,
           18,18]
if args.outType=='xls':
    for i in range(len(colNames)):
        wsw.set_column(i,i,colWidths[i])
    wsw.write_row(0,0,colNames)
elif args.outType=='tsv':
    rFile.write('\t'.join(colNames)+'\n')
    newRows=[]
rowNum=1
for amplName in primerSetOrder:
    try:
        primerSet=primerSets[amplName]
    except KeyError:
        print('ERROR (2)! Unknown amplicon name:')
        print(amplName)
        print('Available amplicon names in the primerSets:')
        print(primerSets.keys())
        print('Type of the searched amplicon name:',type(amplName))
        print('Type of the available amplicon names:',type(list(primerSets.keys())[0]))
        exit(2)
    primers=primerSet.split('_')
    if len(primerInfo[amplName])>4:
        try:
            amplLen=int(primerInfo[amplName][6])-int(primerInfo[amplName][5])+1
        except IndexError:
            print('ERROR! The length of row is less than expected!')
            print(primerInfo[amplName])
            exit(1)
        blockStart=int(primerInfo[amplName][5])+len(primers[0])
        blockEnd=int(primerInfo[amplName][6])-len(primers[1])
    else:
        amplLen=''
        blockStart=''
        blockEnd=''
    lens=[]
    tms=[]
    gcs=[]
    gc_ends=[]
    hairpins=[]
    dimers=[]
    for primer in primers:
        lens.append(len(primer))
        tms.append(round(primer3.calcTm(primer.upper(),
                                        mv_conc=args.mvConc,
                                        dv_conc=args.dvConc,
                                        dntp_conc=0.8,
                                        dna_conc=300,
                                        tm_method='santalucia',
                                        salt_corrections_method='santalucia')
                         ,0))
        gcs.append(round(GC(primer),1))
        gc_ends.append(round(GC(primer[-5:]),1))
        hairpins.append(checkHairpinEnd(primer))
        dimers.append(checkDimerEnd(primer))
    try:
        dimers.append(checkDimerEnd(primers[0],primers[1]))
    except IndexError:
        print('ERROR!',primers)
        print(primerSets)
        exit(1)
    if len(primerInfo[amplName])>4:
        if args.wholeGenomeRef:
            ampliconGC=getAmpliconGC(args.wholeGenomeRef,
                                     primerInfo[amplName][4],
                                     blockStart,
                                     blockEnd)
            ampliconAG=getAmpliconAG(args.wholeGenomeRef,
                                     primerInfo[amplName][4],
                                     int(primerInfo[amplName][5]),
                                     int(primerInfo[amplName][6]))
            ampliconPalindromes=[]
            for mer in range(2,9):
                palindromes=getAmpliconPalindromes(args.wholeGenomeRef,
                                                   primerInfo[amplName][4],
                                                   blockStart,
                                                   blockEnd,
                                                   mer)
                ampliconPalindromes.append(palindromes)
        else:
            ampliconGC=''
            ampliconAG=''
            ampliconPalindromes=['','','','','','','']
        newRow=[*primerInfo[amplName],
                amplLen,blockStart,blockEnd,
                *lens,*tms,*gcs,*gc_ends,
                ampliconGC,ampliconAG,
                *ampliconPalindromes,
                *hairpins,*dimers]
        if args.outType=='xls':
            wsw.write_row(rowNum,0,newRow)
        elif args.outType=='tsv':
            newRow=list(map(str,newRow))
            newRows.append(newRow)
    else:
        newRow=[*primerInfo[amplName],
                '','','',amplLen,blockStart,blockEnd,
                *lens,*tms,*gcs,*gc_ends,
                '','','','','','','','','',
                *hairpins,*dimers]
        if args.outType=='xls':
            wsw.write_row(rowNum,0,newRow)
        elif args.outType=='tsv':
            newRow=list(map(str,newRow))
            newRows.append(newRow)
    rowNum+=1
if args.wholeGenomeRef:
    if len(list(primerInfo.values())[0])>4:
        primerNonSpecRegions=checkPrimersNonTargets(inputFileBase,
                                                    primerSets,
                                                    args.wholeGenomeRef,
                                                    args.substNum,
                                                    args.blastNum,
                                                    args.threads,
                                                    args.maxNonTargetLen,
                                                    args.maxPrimerNonTargets,
                                                    primerInfo)
    else:
        primerNonSpecRegions=checkPrimersNonTargets(inputFileBase,
                                                    primerSets,
                                                    args.wholeGenomeRef,
                                                    args.substNum,
                                                    args.blastNum,
                                                    args.threads,
                                                    args.maxNonTargetLen,
                                                    args.maxPrimerNonTargets,
                                                    None)
    # Now we go through all primer pairs and check them for
    # nonspecific amplicons within one pair of primers
    print(' Searching for nonspecific amplicons that are formed by '
          'input primer pairs...')
    logger.info(' Searching for nonspecific amplicons that are formed by '
                'input primer pairs...')
    wholeWork=len(primerSetOrder)
    for i,amplName in enumerate(primerSetOrder):
        primerSet=primerSets[amplName]
        primerSeq1,primerSeq2=primerSet.split('_')
        primerSpecificity=getNonTargetAmplicons(primerNonSpecRegions,
                                                primerSeq1,
                                                primerSeq2,
                                                args.maxNonTargetLen,
                                                args.threads)
        newRow=[*primerSpecificity[:-1],
                *primerSpecificity[-1]]
        if args.outType=='xls':
            wsw.write_row(i+1,32,newRow)
        elif args.outType=='tsv':
            newRow=newRows[i]+list(map(str,newRow))
            rFile.write('\t'.join(newRow)+'\n')
elif args.outType=='tsv':
    for i in range(len(newRows)):
        rFile.write('\t'.join(newRows[i])+'\n')
if args.outType=='xls':
    wbw.close()
elif args.outType=='tsv':
    rFile.close()
