# This script checks if two primers can form  non-target amplicons

# As input it takes:
# - dictionary (primerNonTargets):
#   {primer_seq:{chromosome:{strand:[position,length]}}} -
#   for primers with mapping
#   {primer_seq:None} - for primers without mappings
# - sequence of the 1st primer that is in the dictionary (primerSeq1)
# - sequence of the 2nd primer that is in the dictionary (primerSeq2)
# - maximal allowed length of non-target amplicon (maxNonTargetLen)

# As output it gives:
# - total number of non-targets for primer1 (totalNonTargetsForPrimer1)
# - total number of non-targets for primer1 (totalNonTargetsForPrimer2)
# - total number of non-target amplicons formed by the primers (nonTargetAmpls)
# - list of all non-target amplicons formed (nonTargetAmplsInRow):
#   [chromosome1,start1,end1,chromosome2,start2,end2,chromosome3,start3,end3...]

import logging
from multiprocessing.pool import ThreadPool
# Own scripts
from work_percentage_process import showPercWork

logger=logging.getLogger(__name__)

def getNonTargetAmplicons(primerNonTargets,
                          primerSeq1,
                          primerSeq2,
                          maxNonTargetLen=200,
                          threads=2):
    if primerSeq1 not in primerNonTargets.keys():
        print('ERROR! The following primer sequence is absent '
              'in the dictionary:')
        logger.error('ERROR! The following primer sequence is absent '
                     'in the dictionary:')
        print(primerNonTargets.keys())
        logger.error(primerSeq1)
        print(primerSeq1)
        exit(1)
    if primerSeq2 not in primerNonTargets.keys():
        print('ERROR! The following primer sequence is absent '
              'in the dictionary:')
        logger.error('ERROR! The following primer sequence is absent '
                     'in the dictionary:')
        print(primerNonTargets.keys())
        logger.error(primerSeq2)
        print(primerSeq2)
        exit(1)
    # Stores number of non-target amplicons formed by these two primers
    nonTargetAmpls=0
    # Stores non-target amplicons formed by these two primers
    nonTargetAmplsInRow=[]
    if (type(primerNonTargets[primerSeq1])==int and
        type(primerNonTargets[primerSeq2])==int):
        totalNonTargetsForPrimer1=primerNonTargets[primerSeq1]
        totalNonTargetsForPrimer2=primerNonTargets[primerSeq2]
        return([totalNonTargetsForPrimer1,
                totalNonTargetsForPrimer2,
                'Unknown',''])
    elif type(primerNonTargets[primerSeq1])==int:
        totalNonTargetsForPrimer1=primerNonTargets[primerSeq1]
        totalNonTargetsForPrimer2=0
    elif type(primerNonTargets[primerSeq2])==int:
        totalNonTargetsForPrimer2=primerNonTargets[primerSeq2]
        totalNonTargetsForPrimer1=0
    else:
        totalNonTargetsForPrimer1=0
        totalNonTargetsForPrimer2=0
    if type(primerNonTargets[primerSeq1])==dict:
        for region in primerNonTargets[primerSeq1].values():
            if 1 in region.keys():
                totalNonTargetsForPrimer1+=len(region[1])
            if -1 in region.keys():
                totalNonTargetsForPrimer1+=len(region[-1])
    if type(primerNonTargets[primerSeq2])==dict:
        for region in primerNonTargets[primerSeq2].values():
            if 1 in region.keys():
                totalNonTargetsForPrimer2+=len(region[1])
            if -1 in region.keys():
                totalNonTargetsForPrimer2+=len(region[-1])
    if (type(primerNonTargets[primerSeq1])==dict and
        type(primerNonTargets[primerSeq2])==dict):
        # Process non-target regions in several threads
        p=ThreadPool(threads)
        # Stores results of processing regions in many threads
        results=[]
##        print('Allocating work between threads...')
##        wholeWork=len(primerNonTargets[primerSeq1])
##        showPercWork(0,wholeWork)
        # We go through only one of dict because
        # we are going to search for nontarget amplicons only 
        ## on chromosomes that are presented in both dicts
        for i,(chrom,
               chrRegions1) in enumerate(primerNonTargets[primerSeq1].items()):
            if chrom not in primerNonTargets[primerSeq2].keys():
                continue
            chrRegions2=primerNonTargets[primerSeq2][chrom]
            results.append(p.apply_async(checkForFormingAmplicon,
                                         (chrRegions1,
                                          chrom,
                                          chrRegions2,
                                          maxNonTargetLen)))
##            showPercWork(i+1,wholeWork)
##        print('\nGetting results...')
        wholeWork=len(results)
        if wholeWork>0:
            showPercWork(0,wholeWork)
            for i,res in enumerate(results):
                res=res.get()
                nonTargetAmpls+=res[0]
                nonTargetAmplsInRow.extend(res[1])
                showPercWork(i+1,wholeWork)
            print()
    return([totalNonTargetsForPrimer1,
            totalNonTargetsForPrimer2,
            nonTargetAmpls,
            nonTargetAmplsInRow])

# This function compares one region to many other regions
# and says if some of them can form amplicon

# As input it takes:
# - coordinates of one region (reg1):
#   {strand:[[position,length],...]}
# - chromosome number (chrom)
# - list of coordinates of the same chromosome
#   to which we compare 1st region (chrRegions2). They have the same format
# - maximal allowed length of non-target amplicon (maxNonTargetLen)

# As output it gives:
# - total number of non-target amplicons formed by the primers (nonTargetAmpls)
# - list of all non-target amplicons formed (nonTargetAmplsInRow):
#   [chromosome1,start1,end1,chromosome2,start2,end2,chromosome3,start3,end3...]
def checkForFormingAmplicon(chrRegions1,chrom,chrRegions2,maxNonTargetLen):
    nonTargetAmpls=0
    nonTargetAmplsInRow=[]
    # Initially, we search for amplicons that are formed by 1st primer on plus
    # strand and 2nd primer on minus strand
    if (1 in chrRegions1.keys() and
        -1 in chrRegions2.keys()):
        bothChrRegions=[]
        for reg1 in chrRegions1[1]:
            bothChrRegions.append([1,reg1])
        for reg2 in chrRegions2[-1]:
            bothChrRegions.append([-1,reg2])
        bothChrRegions=sorted(bothChrRegions,
                              key=lambda item:item[1][0])
        for i,reg1 in enumerate(bothChrRegions):
            if reg1[0]<0:
                continue
            for reg2 in bothChrRegions[i+1:]:
                if reg2[0]>0:
                    continue
                if reg2[1][0]+reg2[1][1]-reg1[1][0]>maxNonTargetLen:
                    break
                if 0<reg2[1][0]+reg2[1][1]-reg1[1][0]<=maxNonTargetLen:
                    nonTargetAmpls+=1
                    nonTargetAmplsInRow.extend([chrom,
                                                reg1[1][0],
                                                reg2[1][0]+reg2[1][1]-1])
    # Then, we search for amplicons that are formed by 1st primer on minus
    # strand and 2nd primer on plus strand
    if (-1 in chrRegions1.keys() and
        1 in chrRegions2.keys()):
        bothChrRegions=[]
        for reg1 in chrRegions1[-1]:
            bothChrRegions.append([-1,reg1])
        for reg2 in chrRegions2[1]:
            bothChrRegions.append([1,reg2])
        bothChrRegions=sorted(bothChrRegions,
                              key=lambda item:item[1][0])
        # Here reg1 and reg2 - regions on plus and minus strands, respectively
        for i,reg1 in enumerate(bothChrRegions):
            if reg1[0]<0:
                continue
            for reg2 in bothChrRegions[i+1:]:
                if reg2[0]>0:
                    continue
                if reg2[1][0]+reg2[1][1]-reg1[1][0]>maxNonTargetLen:
                    break
                if 0<reg2[1][0]+reg2[1][1]-reg1[1][0]<=maxNonTargetLen:
                    nonTargetAmpls+=1
                    nonTargetAmplsInRow.extend([chrom,
                                                reg1[1][0],
                                                reg2[1][0]+reg2[1][1]-1])
    return(nonTargetAmpls,nonTargetAmplsInRow)
