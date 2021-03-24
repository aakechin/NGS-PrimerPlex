# This function reads BAM-file of mapping primer sequence onto reference genome
# and returns dictionary that contains all regions to which each primer sequence
# was mapped

# As input it gets:
# - AlignedSegment from pysam (read)
# - maximum number of nonspecific regions that user is ready to save (maxPrimerNonspec)
# - reference genome opened with pysam (refFa)
# - information about primers for storing target region (primersInfo):
#   {amplicon_name:[amplicon_number,left_primer,right_primer,
#                   amplicon_name,amplicon_start,amplicon_end]}

# As output it gives:
# - sequence of primer that was mapped onto the reference genome
# - list of unspecific regions to which it was mapped:
#   [[chrom,strand,position,regionLen]]

import re
import logging
from Bio import Seq
from Bio import pairwise2

logger=logging.getLogger(__name__)

def revComplement(nuc):
    try:
        return(str(Seq.Seq(nuc).reverse_complement()))
    except ValueError as e:
        if 'Mixed RNA/DNA found' in str(e):
            print('ERROR (3): Mixed RNA/DNA found:')
            print(nuc)
            exit(3)
        else:
            print('ERROR: Unknown ValueError!')
            exit(0)

def readToPrimerNonTargets(read,maxPrimerNonspec,refFa,primersInfo=None):
    indelsPat=re.compile('(\d+)([ID])')
    matchPat=re.compile('(\d+)M')
    # Extract strand of the main match in genome
    if read.flag==0:
        strand=1
    elif read.flag==16 or read.flag==20:
        strand=-1
    elif read.flag==4:
        return(read.qname,[])
    else:
        print('ERROR! Unknown value of FLAG:',read.flag)
        print(read)
        exit(1)
    if primersInfo:
        primersName=read.qname
        for primerPairNum,primerPairInfo in primersInfo.items():
            if primersName in primerPairInfo[1:3]:
                primerNumInPair=primerPairInfo[1:3].index(primersName)
                infoStrand=int((-1)**primerNumInPair)
                try:
##                    pos=int(primersInfo[primerPairNum][5+primerNumInPair*4])+primerNumInPair
                    pos=int(primersInfo[primerPairNum][5+primerNumInPair])
                except IndexError:
                    print('ERROR (2): Incorrect index:')
                    print(primersInfo)
                    print(primersName)
                    print(primerPairNum)
                    print(primerNumInPair)
                    exit(2)
                try:
                    targetRegion=[primersInfo[primerPairNum][4],
                                  infoStrand*pos]
                except TypeError:
                    print('ERROR!',primersInfo[primerPairNum][4],primersInfo[primerPair][0],primerNumInPair)
                    exit(15)
                break
    else:
        targetRegion=[]
    if strand==-1:
        mainMapping=[read.reference_name,strand*(read.pos+len(read.qname))]
    else:
        mainMapping=[read.reference_name,strand*(read.pos+1)]
    if (mainMapping!=targetRegion):
        ## The next string is necessary for human genomes
        ## to exclude unsorted chromosome fragments from the analysis
##        and '_' not in read.reference_name
        nonSpecRegions=[','.join([read.reference_name,
                                  str(strand*(read.pos+1)),
                                  read.cigarstring,
                                  str(read.get_tag('NM'))])]
    else:
        nonSpecRegions=[]
    if read.has_tag('XA'):
        # XA tag of read ends with ; so the last element is empty
        nonSpecRegions.extend(read.get_tag('XA').split(';')[:-1])
    qname=read.qname
    if len(nonSpecRegions)>maxPrimerNonspec:
        return(qname,nonSpecRegions)
    primerNonSpecRegions=[]
    for region in nonSpecRegions: 
        # We go through all these regions and check
        # that 3'-nucleotide matches primer's 3'-end
        chrom,pos,cigar,subst=region.split(',')
        # Determine length of sequence of interest by parsing CIGAR.
        # The length depends only on the deletions
        # but not mismatches nor insertions into reference genome
        # So we count number of deletions
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
            pos2=-(int(pos)+regionLen)
        else:
            pos2=int(pos)
        if ([chrom,int(pos)]==targetRegion or
            (len(targetRegion)>0 and
             chrom==targetRegion[0] and
             abs(targetRegion[1]-pos2)<=len(read.seq))):
##            print([chrom,int(pos)])
            continue
        attempts=0
        seq=None
        while(seq is None):
            try:
                seq=refFa.fetch(region=chrom+':'
                                ''+str(abs(int(pos)))+'-'
                                ''+str(abs(int(pos))+regionLen-1))
            except:
                seq=None
                attempts+=1
                if attempts>=10:                    
                    print('ERROR!')
##                    logger.error(str(e))
                    print(refFa.filename)
                    logger.error(refFa.filename)
                    print(chrom,pos,regionLen)
                    exit(1)
        # Determine, which sequence we should take:
        # forward or reverse-complement
        # If primer is on + strand and found region is on opposite
        # or primer is on - strand and found region is on the same
        # we take reverse-complement
        if int(pos)>0:
            regStrand=1
        elif int(pos)<0:
            regStrand=-1
        if (strand>0 and int(pos)<0) or (strand<0 and int(pos)>0):
            try:
                seq=str(Seq.Seq(seq).reverse_complement()).upper()
            # If there is an error like 'Mixed RNA/DNA found' 
            except ValueError as e:
                if 'Mixed RNA/DNA found' in str(e):
                    print('ERROR (1): Mixed RNA/DNA found:')
                    print(seq)
                    exit(1)
                else:
                    print('ERROR: Unknown ValueError!')
                    exit(0)
        else:
            seq=seq.upper()
        # Check if read sequence is the same as read qname
        if read.qname==read.seq:
            primerSeq=read.seq
            regionSeq=seq
        else:
            primerSeq=read.qname
            regionSeq=revComplement(seq)
        if len(regionSeq)==1:
            continue
        # If there is some insertions or deletion in found region
        if 'I' in cigar or 'D' in cigar:
            # We need to align its sequence with primer sequence
            align=pairwise2.align.globalxx(primerSeq,regionSeq)
            # Check that 3'-ends of primers are identical
            ## and one of two nucleotides before 3'-ends are identical, too
            if (align[0][0][-1]==align[0][1][-1] or
                align[0][0][-2]==align[0][1][-2]):
                # Then we consider this region as a non-specific
                # for this primer
                primerNonSpecRegions.append([chrom,
                                             regStrand,
                                             abs(int(pos)),
                                             regionLen])
        else:
            try:
                if (regionSeq[-1]==primerSeq[-1] or
                    regionSeq[-2]==primerSeq[-2]):
                    # Then we consider this region as a non-specific for this primer
                    primerNonSpecRegions.append([chrom,
                                                 regStrand,
                                                 abs(int(pos)),
                                                 regionLen])
            except IndexError:
                print('ERROR (2): incorrect index for sequences:')
                print('regionSeq:',regionSeq)
                print('primerSeq:',primerSeq)
                print('chr:',chrom)
                print('position:',pos)
                print('regionLen:',regionLen)
                print(-1,-2)
                exit(2)
    return(read.qname,primerNonSpecRegions)
