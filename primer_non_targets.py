# This script checks primer pairs for unspecific hybridization in genome

# As input it takes:
# - base for the FASTA-file name with primer sequences (inputFileBase)
# - all primer pairs (primerSets):
#   {amplicon_name:seq1_seq2}
# - path to the indexed reference genome sequence (wholeGenomeRef)
# - number of substitutions which are allowed during mapping (substNum)
# - number of threads of reading BAM-file (threads)
# - maximum length of an unspecific amplicon (maxNonSpecLen)
# - maximum number of nonspecific regions to which primer can be mapped
#   and we will still search for nonspecific amplicons (maxPrimerNonspec)
# - information about primers for storing target region (primersInfo):
#   {amplicon_name:[amplicon_number,left_primer,right_primer,
#                   amplicon_name,amplicon_start,amplicon_end]}

# As output it gives:
# - dictionary (primerNonTargets):
#   {primer_seq:{chromosome:{strand:[position,length]}}} -
#   for primers with mapping
#   {primer_seq:None} - for primers without mappings

import logging
import pysam
import subprocess as sp
from multiprocessing.pool import ThreadPool
# Own scripts
from primer_bam_to_non_targets import readToPrimerNonTargets
from work_percentage_process import showPercWork

logger=logging.getLogger(__name__)

def checkPrimersNonTargets(inputFileBase,primerSets,wholeGenomeRef,
                            substNum=1,bnum=1,threads=2,maxNonSpecLen=100,
                            maxPrimerNonspec=1000,primersInfo=None):
    # Dictionary for storing info about primers specificity by primers
    # for checking specificity within one amplicon    
    primerNonTargets={}
    # Creating fasta-file with all primers' sequences
    seqFile=open(inputFileBase+'_checkPrimers_all_primers_sequences.fa','w')
    for amplName,key in primerSets.items():
        primers=key.split('_')
        for j,primer in enumerate(primers):
            if primer!='' and primer not in primerNonTargets.keys():
                seqFile.write('\n'.join(['>'+primer,primer])+'\n')
                primerNonTargets[primer]=None
            elif primer not in primerNonTargets.keys():
                primerNonTargets[primer]=None
    seqFile.close()
    bwaResultFileName=inputFileBase+'_checkPrimers_all_primers_sequences.bwa'
    for runNum in range(bnum):
        if runNum==1:
            print('\n Running BWA again to search all regions...')
            logger.info(' Running BWA again to search all regions...')
        else:
            print(' Running BWA...')
            logger.info(' Running BWA...')
        out=sp.check_output(' '.join(['bwa','aln',
                                      '-N','-n',str(substNum),
                                      '-t',str(threads),wholeGenomeRef,
                                      seqFile.name,'>',
                                      bwaResultFileName+'.sai']),
                            shell=True,stderr=sp.STDOUT).decode('utf-8')
        out=sp.check_output(' '.join(['bwa','samse',
                                      '-n','100000000',
                                      wholeGenomeRef,
                                      bwaResultFileName+'.sai',
                                      seqFile.name,'>',
                                      bwaResultFileName+'.sam']),
                                      shell=True,stderr=sp.STDOUT).decode('utf-8')
        # Reading BWA output file
        samFile=pysam.AlignmentFile(bwaResultFileName+'.sam')
        # Process SAM-file strings in several threads
        p=ThreadPool(threads)
        print(' Processing SAM-file...')
        logger.info(' Processing SAM-file...')
        results=[]
        refFa=pysam.FastaFile(wholeGenomeRef)
        for read in samFile.fetch():
            results.append(p.apply_async(readToPrimerNonTargets,(read,
                                                                 maxPrimerNonspec,
                                                                 refFa,
                                                                 primersInfo)))
        wholeWork=len(results)
        doneWork=0
        for res in results:
            # each res is:
            # - sequence of primer that was mapped onto the reference genome
            # - list of unspecific regions to which it was mapped:
            #   [[chrom,strand,position,regionLen]]
            res=res.get()
            doneWork+=1
            showPercWork(doneWork,wholeWork)
            if res[1]==None or len(res[1])==0:
                primerNonTargets[res[0]]=None
            else:
                if primerNonTargets[res[0]]==None:
                    primerNonTargets[res[0]]={}
                if len(res[1])>maxPrimerNonspec:
                    primerNonTargets[res[0]]=len(res[1])
                    continue
                # res[1] has format: chrom, strand, pos, length
                for region in res[1]:
                    if region[0] not in primerNonTargets[res[0]].keys():
                        primerNonTargets[res[0]][region[0]]={region[1]:set([tuple(region[2:])])}
                    elif region[1] not in primerNonTargets[res[0]][region[0]].keys():
                        primerNonTargets[res[0]][region[0]][region[1]]=set([tuple(region[2:])])
                    else:
                        primerNonTargets[res[0]][region[0]][region[1]].add(tuple(region[2:]))
        p.close()
        p.join()
        refFa.close()
        totalNumberNonTargets=0
        for primer,nonTargets in sorted(primerNonTargets.items()):
            if type(nonTargets)==int:
                sumNum=nonTargets
            else:
                sumNum=0
                if nonTargets!=None:
                    for chrom,strandRegions in nonTargets.items():
                        for regions in strandRegions.values():
                            sumNum+=len(regions)
            totalNumberNonTargets+=sumNum
    print('\n # Total number of nonspecific regions:',
          totalNumberNonTargets)
    logger.info(' '.join([' # Total number of nonspecific regions:',
                          str(totalNumberNonTargets)]))
    return(primerNonTargets)
