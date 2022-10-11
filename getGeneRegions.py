# This script takes names of genes and numbers of their exons or positions in CDS
# and makes region file for NGS-PrimerPlex program

# Section of importing modules
from Bio import SeqIO
import copy
import argparse
import os
import glob
import re
import gzip
import pysam

# Global variables
global thisDir,nameToNum,numToName
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

# Section of functions
def chrToChr(wholeGenomeRef):
    global nameToNum,numToName
    nameToNum={}
    numToName={}
    try:
        refFa=pysam.FastaFile(wholeGenomeRef)
    except OSError:
        print('ERROR (16)! Reference genome FASTA-file is absent:')
        print(wholeGenomeRef)
        print(glob.glob(os.path.dirname(wholeGenomeRef)))
        exit(16)
    for i,ch in enumerate(refFa.references):
        nameToNum[ch]=i+1
        numToName[i+1]=ch
    return(nameToNum,numToName)

def getChrNum(geneName,refDir):
    global nameToNum,numToName
    if not os.path.exists(refDir+'geneNameToChromosome.csv'):
        print('WARNING! There was no table for converting gene name into chromosome.')
        print('It will be created and this process will take some time...')
        try:
            geneNameToChromosomes=createGeneNameToChromosomeFile(refDir)
        except ValueError as e:
            if os.path.exists(refDir+'geneNameToChromosome.csv'):
                os.remove(refDir+'geneNameToChromosome.csv')
            if 'More than one record found in handle' in str(e):
                print('ERROR (17)! Table for converting gene name into chromosome could not be created ')
                print('because some of GenBank files used contain more than one fragment. Each fragment '
                      '(e.g. chromosome or plasmid) should be written into distinct GenBank-file')
                exit(17)
            else:
                print('ERROR (12)! Table for converting gene name into chromosome could not be created')
                exit(12)
        except:
            if os.path.exists(refDir+'geneNameToChromosome.csv'):
                os.remove(refDir+'geneNameToChromosome.csv')
            print('ERROR (12)! Table for converting gene name into chromosome could not be created')
            exit(12)
        print('Table for converting gene name into chromosomes was created:')
        print(refDir+'geneNameToChromosome.csv')
    else:
        geneNameToChromosomes=readGeneNameToChromosome(refDir+'geneNameToChromosome.csv')
    if geneName in geneNameToChromosomes.keys():
        if len(geneNameToChromosomes[geneName])==1:
            return(next(iter(geneNameToChromosomes[geneName])))
        else:
            print('ERROR (13)! The following gene has several locations in genome:')
            print(geneName)
            print('Chromosomes: '+', '.join(geneNameToChromosomes[geneName]))
            print('Use another synonym for these gene (find it in the GenBank-file)')
            exit(13)
    else:
        print('ERROR (14)! The following gene was not found in the reference genome:')
        print(geneName)
        print('Reference genome location: '+refDir)
        exit(14)

def readGeneNameToChromosome(fileName):
    geneNameToChromosomes={}
    with open(fileName) as file:
        for string in file:
            cols=string.replace('\n','').split('\t')
            geneNameToChromosomes[cols[0]]=set()
            for chrom in cols[1:]:
                if chrom!='':
                    geneNameToChromosomes[cols[0]].add(chrom)
    return(geneNameToChromosomes)

def createGeneNameToChromosomeFile(refDir):
    ds=glob.glob(refDir+'*.gb')
    if len(ds)==0:
        ds=glob.glob(refDir+'*.gb.gz')
        if len(ds)==0:
            print('ERROR (10)! In the defined reference directory '
                  'there is no GenBank files that starts from chr')
            print(refDir)
            exit(10)
    geneNameToChromosomes={}
    for d in sorted(ds):
        print(os.path.basename(d))
        chromNumOrName=os.path.splitext(os.path.basename(d))[0]
        try:
            chromNum=int(chromNumOrName)
            if int(chromNum) in numToName.keys():
                chrom=numToName[int(chromNum)]
            else:
                print('ERROR (15)! The following GenBank-file has incorrect format of name:')
                print(os.path.basename(d))
                print('It should have name in the format of the reference genome FASTA-file, e.g. chr1.gb')
                print('Or all files can be named with numbers of chromosomes as they are written in the reference genome FASTA-file.')
                print('In this case rename e.g. for hg19 they will be:')
                print('chrM - 1.gb')
                print('chr1 - 2.gb')
                print('chr2 - 3.gb')
                print('. . .')
                print('chrX - 24.gb')
                print('chrY - 25.gb')
                exit(15)
        except:
            chromName=chromNumOrName
            if chromName in nameToNum.keys():
                chrom=chromName
            else:
                print('ERROR (9)! The following GenBank-file has incorrect format of name:')
                print(os.path.basename(d))
                print('It should have name in the format of the reference genome FASTA-file, e.g. chr1.gb')
                print('Or all files can be named with numbers of chromosomes as they are written in the reference genome FASTA-file.')
                print('In this case rename e.g. for hg19 they will be:')
                print('chrM - 1.gb')
                print('chr1 - 2.gb')
                print('chr2 - 3.gb')
                print('. . .')
                print('chrX - 24.gb')
                print('chrY - 25.gb')
                exit(9)
        if d[-3:]=='.gz':
            data=SeqIO.read(gzip.open(d,'rt'),'genbank')
        else:
            data=SeqIO.read(d,'genbank')
        for f in data.features:
            if f.type!='gene':
                continue
            if 'gene_synonym' in f.qualifiers.keys():
                synonyms=f.qualifiers['gene_synonym'][0].split('; ')
            else:
                synonyms=[]
            if 'gene' in f.qualifiers.keys():
                geneName=f.qualifiers['gene'][0]
                if geneName not in geneNameToChromosomes.keys():
                    geneNameToChromosomes[geneName]=set()
                geneNameToChromosomes[geneName].add(chrom)
                for synonym in synonyms:
                    if synonym not in geneNameToChromosomes.keys():
                        geneNameToChromosomes[synonym]=set()
                    geneNameToChromosomes[synonym].add(chrom)
    writeGeneNameToChromosomeFile(geneNameToChromosomes,refDir)
    return(geneNameToChromosomes)

def writeGeneNameToChromosomeFile(geneNameToChromosomes,refDir):
    with open(refDir+'geneNameToChromosome.csv','w') as file:
        for geneName,chroms in geneNameToChromosomes.items():
            file.write('\t'.join([geneName,*chroms])+'\n')

def writeRegions(resultFile,targetGene,pro_mRNA,cds=None,codons=[],exons=[],nucs=[]):
    # Convert list of exons coordinates to dictionary with exon numbers as keys and coordinates as values
    mRNA={}
    for exonNum,part in enumerate(pro_mRNA):
        mRNA[exonNum]=[part.nofuzzy_start,part.nofuzzy_end,part.strand]
##    mRNA=dict(zip(range(len(mRNA)),mRNA))
    # If we don't include noncoding regions, we need to change
    if cds:
        coding_mRNA=mRNA.copy()
        # We need to remove noncoding exons or exon parts
        if cds[0].strand>0:
            codingStart=cds[0].nofuzzy_start+1
            codingEnd=cds[-1].nofuzzy_end
            for exonNum in sorted(mRNA.keys()):
                # If whole exon is located before start of CDS, we remove it
                if mRNA[exonNum][1]<codingStart:
                    coding_mRNA.pop(exonNum)
                elif mRNA[exonNum][0]+1<codingStart:
                    coding_mRNA[exonNum][0]=codingStart-1
                    break
                else: break
            for exonNum in sorted(mRNA.keys())[::-1]:
                # If whole exon is located after end of CDS, we remove it
                if mRNA[exonNum][0]+1>codingEnd:
                    coding_mRNA.pop(exonNum)
                elif mRNA[exonNum][1]>codingEnd:
                    coding_mRNA[exonNum][1]=codingEnd
                    break
                else: break
        elif cds[0].strand<0:
            codingStart=cds[0].nofuzzy_end
            codingEnd=cds[-1].nofuzzy_start+1
            for exonNum in sorted(mRNA.keys()):
                # If whole exon is located after start of CDS, we remove it
                if mRNA[exonNum][0]+1>codingStart:
                    coding_mRNA.pop(exonNum)
                elif mRNA[exonNum][1]>codingStart:
                    coding_mRNA[exonNum][1]=codingStart
                    break
                else: break
            for exonNum in sorted(mRNA.keys())[::-1]:
                # If whole exon is located before start of CDS, we remove it
                if mRNA[exonNum][1]<codingEnd:
                    coding_mRNA.pop(exonNum)
                elif mRNA[exonNum][0]+1<codingEnd:
                    coding_mRNA[exonNum][0]=codingEnd-1
                    break
                else: break
        mRNA=coding_mRNA.copy()
    if (len(codons)==0 and
        len(exons)==0):
        for exonNum,subf in mRNA.items():
            if exonNum+1!=1 and exonNum+1!=len(mRNA):
                resultFile.write('\t'.join([chrs[targetGene],
                                            str(subf[0]+1-intronSize),
                                            str(subf[1]+intronSize),
                                            '_'.join([targetGene,
                                                          'ex'+str(exonNum+1)]),
                                            '1','B','NotW'])+'\n')
            elif (exonNum+1==1 and subf[2]>0) or (exonNum+1==len(mRNA) and subf[2]<0):
                resultFile.write('\t'.join([chrs[targetGene],
                                            str(subf[0]+1),
                                            str(subf[1]+intronSize),
                                            '_'.join([targetGene,
                                                          'ex'+str(exonNum+1)]),
                                            '1','B','NotW'])+'\n')
            elif (exonNum+1==1 and subf[2]<0) or (exonNum+1==len(mRNA) and subf[2]>0):
                resultFile.write('\t'.join([chrs[targetGene],
                                            str(subf[0]+1-intronSize),
                                            str(subf[1]),
                                            '_'.join([targetGene,
                                                          'ex'+str(exonNum+1)]),
                                            '1','B','NotW'])+'\n')
            else:
                print('ERROR (1)! Unknown variant of exon number and gene strand!')
                print(mRNA)
                print(exonNum+1,subf[2])
                exit(1)
    else:
        if len(codons)>0:
            exLen=0
            exLens={}
            # If the 1st coding exon is no the 1st exon
            for exonNum in range(min(mRNA.keys())):
                exLens[exonNum]=0
            for exonNum,subf in sorted(mRNA.items()):
                exLen+=subf[1]-subf[0]
                exLens[exonNum]=exLen
            for regNum,reg in enumerate(nucs):
                for z,r in enumerate(reg):
                    exFound=False
                    for j,exLen in exLens.items():
                        if r<=exLen:
                            exFound=True
                            # If it is not the first exon
                            if j!=0:
                                # if it is start of necessary region
                                if z==0:
                                    fromStart1=r-exLens[j-1]-1
                                    exNum1=j
                                else:
                                    fromStart2=r-exLens[j-1]-1
                                    exNum2=j
                            else:
                                if z==0:
                                    fromStart1=r
                                    exNum1=j
                                else:
                                    fromStart2=r
                                    exNum2=j
                            break
                if not exFound:
                    print('ERROR (2)! Exon for region '+str(reg)+' of gene '+targetGene+' was not determined correctly')
                    exit(2)
                if codons[regNum][0]==codons[regNum][1]:
                    regStr=str(codons[regNum][0])
                else:
                    regStr='_'.join(list(map(str,codons[regNum])))
                # If exons for the start and for the end of region are equal
                if exNum1==exNum2:
                    if list(mRNA.values())[0][2]==-1:
                        try:
                            resultFile.write('\t'.join([chrs[targetGene],
                                                        str(mRNA[exNum2][1]-fromStart2),
                                                        str(mRNA[exNum1][1]-fromStart1),
                                                        '_'.join([targetGene,'p'+regStr]),
                                                        '1','B','NotW'])+'\n')
                        except KeyError:
                            print('ERROR (3)!',exNum1,exNum2)
                            print(mRNA)
                            exit(3)
                    else:
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum1][0]+1+fromStart1),
                                                    str(mRNA[exNum2][0]+1+fromStart2),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                elif exNum2-exNum1==1:
                    if list(mRNA.values())[0][2]==-1:
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum1][0]+1),
                                                    str(mRNA[exNum1][1]-fromStart1),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum2][1]-fromStart2),
                                                    str(mRNA[exNum2][1]),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                    else:
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum1][0]+1+fromStart1),
                                                    str(mRNA[exNum1][1]),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum2][0]+1),
                                                    str(mRNA[exNum2][0]+1+fromStart2),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                else:
                    if list(mRNA.values())[0][2]==-1:
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum1][0]+1),
                                                    str(mRNA[exNum1][1]-fromStart1),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                        for j in range(exNum1+1,exNum2):
                            resultFile.write('\t'.join([chrs[targetGene],
                                                        str(mRNA[j][0]+1),
                                                        str(mRNA[j][1]),
                                                        '_'.join([targetGene,'p'+regStr]),
                                                        '1','B','NotW'])+'\n')
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum2][1]-fromStart2),
                                                    str(mRNA[exNum2][1]),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                    else:
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum1][0]+1+fromStart1),
                                                    str(mRNA[exNum1][1]),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
                        for j in range(exNum1+1,exNum2):
                            resultFile.write('\t'.join([chrs[targetGene],
                                                        str(mRNA[j][0]+1),
                                                        str(mRNA[j][1]),
                                                        '_'.join([targetGene,'p'+regStr]),
                                                        '1','B','NotW'])+'\n')
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(mRNA[exNum2][0]+1),
                                                    str(mRNA[exNum2][0]+1+fromStart2),
                                                    '_'.join([targetGene,'p'+regStr]),
                                                    '1','B','NotW'])+'\n')
        if len(exons)>0:
            for exRange in exons:
                for exonNum in range(exRange[0]-1,exRange[1]):
                    if exonNum not in mRNA.keys():
                        print('ERROR (4)! Exon '+str(exonNum)+' includes noncoding sequences or does not exist.')
                        print('In the first case, use parameter -noncoding. In the second one, correct the mistake')
                        exit(4)
                    subf=mRNA[exonNum]
##                    resultFile.write('\t'.join([chrs[targetGene],
##                                                str(subf[0]+1),
##                                                str(subf[1]),
##                                                '_'.join([targetGene,
##                                                          'ex'+str(exonNum+1)]),
##                                                '1','B','NotW'])+'\n')
                    if exonNum+1!=1 and exonNum+1!=len(mRNA):
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(subf[0]+1-intronSize),
                                                    str(subf[1]+intronSize),
                                                    '_'.join([targetGene,
                                                          'ex'+str(exonNum+1)]),
                                                    '1','B','NotW'])+'\n')
                    elif (exonNum+1==1 and subf[2]>0) or (exonNum+1==len(mRNA) and subf[2]<0):
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(subf[0]+1),
                                                    str(subf[1]+intronSize),
                                                    '_'.join([targetGene,
                                                          'ex'+str(exonNum+1)]),
                                                    '1','B','NotW'])+'\n')
                    elif (exonNum+1==1 and subf[2]<0) or (exonNum+1==len(mRNA) and subf[2]>0):
                        resultFile.write('\t'.join([chrs[targetGene],
                                                    str(subf[0]+1-intronSize),
                                                    str(subf[1]),
                                                    '_'.join([targetGene,
                                                          'ex'+str(exonNum+1)]),
                                                    '1','B','NotW'])+'\n')
                    else:
                        print('ERROR (11)! Unknown variant of exon number and gene strand!')
                        print(mRNA)
                        print(exonNum+1,subf[2])
                        exit(1)

# Section of reading arguments
parser=argparse.ArgumentParser(description='This script takes names of genes '
                               'and numbers of their exons or positions in CDS '
                               'and makes region file for hi-plex program')
parser.add_argument('--geneListFile','-glf',type=str,
                    help='file with list of genes. Format is: GENE EXONS CODONS',required=True)
parser.add_argument('--refDir','-ref',type=str,
                    help='directory with reference files',required=True)
parser.add_argument('--reference-genome','-wgref',
                 dest='wholeGenomeRef',type=str,
                 help='file with INDEXED whole-genome reference sequence',
                 required=True)
parser.add_argument('--resultFile','-rf',type=str,
                    help='file for results',required=True)
parser.add_argument('--intron-nucleotides','-intron',type=int,dest='intronSize',
                    help='number of nucleotides from intron to take. Default: 2',default=2)
parser.add_argument('--include-noncoding','-noncoding',dest='includeNonCoding',action='store_true',
                    help="use this parameter, if you want to include 5'- and 3'-non-coding regions of mRNA")
args=parser.parse_args()
chrToChr(args.wholeGenomeRef)
refDir=args.refDir
intronSize=args.intronSize

if refDir[-1]!=os.path.sep:
    refDir+=os.path.sep                

# Read input file
print('Reading input file and getting chromosome numbers...')
targetGenes=[]
chrs={}
codons=[]
exons=[]
geneListFile=open(args.geneListFile)
for string in geneListFile:
    if 'Gene\tExons' in string:
        continue
    cols=string[:-1].split('\t')
    print(cols[0])
    targetGenes.append(cols[0])
    exons.append([])
    codons.append([])
    chrs[cols[0]]=getChrNum(cols[0],refDir)
    if chrs[cols[0]]==None:
        print('ERROR (5)! Gene '+cols[0]+' was not found in the Gene database of NCBI!')
        exit(5)
    if not (len(cols)<2 or cols[1]==''):
        ex=cols[1].split(',')
        for e in ex:
            if '-' in e:
                ep=e.split('-')
                if [int(ep[0]),int(ep[1])] not in exons[-1]:
                    exons[-1].append([int(ep[0]),int(ep[1])])
            else:
                if [int(e),int(e)] not in exons[-1]:
                    exons[-1].append([int(e),int(e)])
    if not (len(cols)<3 or cols[2]==''):
        ex=cols[2].split(',')
        for e in ex:
            if '-' in e:
                ep=e.split('-')
                if [int(ep[0]),int(ep[1])] not in codons[-1]:
                    codons[-1].append([int(ep[0]),int(ep[1])])
            else:
                if [int(e),int(e)] not in codons[-1]:
                    codons[-1].append([int(e),int(e)])
geneListFile.close()
print('Done')

# Library nucs contains ranges of nucleotides for which we need to line up primers
# we get it from codons list
nucs=copy.deepcopy(codons)
for key,value in enumerate(nucs):
    for j,l in enumerate(value):
        nucs[key][j][0]=l[0]*3-2
        nucs[key][j][1]=l[1]*3
resultFile=open(args.resultFile,'w')
transcriptPat=re.compile('transcript variant ([A-z\d]+)')
# For some genes word "isoform" is repeated in a GB-file
isoformPat=re.compile('(?:isoform )+([X\da-z]+)')
# Dictionary that stores genes and their main transcript accession numbers
mainTrascripts={'NRG1':'NM_013957.'}
mainProteins={'NRG1':'NP_039251.'}
## TO DO:
# Maybe add this values to input table for such complex genes?
print('Getting genome regions for the input genes...')
for geneNum,t in enumerate(targetGenes):
    print(t)
    # Here, chrs[t] will have chromosome name,
    # but GB-files can be named by chromosome numbers
    # Names
    if os.path.exists(refDir+str(chrs[t])+'.gb'):
        data=SeqIO.read(refDir+str(chrs[t])+'.gb','genbank')
    elif os.path.exists(refDir+str(chrs[t])+'.gb.gz'):
        data=SeqIO.read(gzip.open(refDir+str(chrs[t])+'.gb.gz','rt'),'genbank')
    # Numbers
    elif os.path.exists(refDir+str(nameToNum[chrs[t]])+'.gb'):
        data=SeqIO.read(refDir+str(nameToNum[chrs[t]])+'.gb','genbank')
    elif os.path.exists(refDir+str(nameToNum[chrs[t]])+'.gb.gz'):
        data=SeqIO.read(gzip.open(refDir+str(nameToNum[chrs[t]])+'.gb.gz','rt'),'genbank')
    else:
        print('ERROR (6)! GenBank file for chromosome '+str(chrs[t])+' is absent in the reference directory '+refDir+'!')
        print(refDir+str(chrs[t])+'.gb', refDir+str(chrs[t])+'.gb.gz', str(nameToNum[chrs[t]])+'.gb','or',str(nameToNum[chrs[t]])+'.gb.gz')
        exit(6)
    geneFound=False
    # several_mRNA_isoforms contains all saved mRNA isoforms for current gene
    ## And then we take one with the lowest number
    several_mRNA_isoformsNums={}
    several_mRNA_isoformsWords={}
    mRNA_exons=[]
    for f in data.features:
        if 'gene_synonym' in f.qualifiers.keys():
            synonyms=f.qualifiers['gene_synonym'][0].split('; ')
        if len(mRNA_exons)==0 and ((('gene' in f.qualifiers.keys() and f.qualifiers['gene'][0]==t) or ('gene_synonym' in f.qualifiers.keys() and t in synonyms)) and
            f.type=='mRNA' and ('product' in f.qualifiers.keys())):
            if (t in mainTrascripts.keys() and
                'transcript_id' in f.qualifiers.keys() and
                mainTrascripts[t] in f.qualifiers['transcript_id'][0]):
                mRNA_exons=f.location.parts
                if args.includeNonCoding and t not in codons.keys():
                    geneFound=True
                    writeRegions(resultFile,t,mRNA_exons,None,codons[geneNum],exons[geneNum],nucs[geneNum])
                    break
            elif (len(transcriptPat.findall(f.qualifiers['product'][0]))==0 or
                transcriptPat.findall(f.qualifiers['product'][0])[0]=='1' or
                transcriptPat.findall(f.qualifiers['product'][0])[0]=='a'):
                mRNA_exons=f.location.parts
                if args.includeNonCoding and t not in codons.keys():
                    geneFound=True
                    writeRegions(resultFile,t,mRNA_exons,None,codons[geneNum],exons[geneNum],nucs[geneNum])
                    break
            elif len(transcriptPat.findall(f.qualifiers['product'][0]))>0:
                try:
                    key=int(transcriptPat.findall(f.qualifiers['product'][0])[0])
                    several_mRNA_isoformsNums[key]=f.location.parts
                except ValueError:
                    key=transcriptPat.findall(f.qualifiers['product'][0])[0]
                    several_mRNA_isoformsWords[key]=f.location.parts
        if ((('gene' in f.qualifiers.keys() and
              f.qualifiers['gene'][0]==t) or
             ('gene_synonym' in f.qualifiers.keys() and
              t in synonyms)) and
            f.type=='CDS'):
            if 'product' in f.qualifiers.keys():
                isoformMatch=isoformPat.findall(f.qualifiers['product'][0])
            else:
                isoformMatch=[]
            if 'note' in f.qualifiers.keys():
                transcriptMatch=transcriptPat.findall(f.qualifiers['note'][0])
            else:
                transcriptMatch=[]
            if ((t in mainProteins.keys() and
                 'protein_id' in f.qualifiers.keys() and
                 mainProteins[t] in f.qualifiers['protein_id'][0]) or
                ('product' in f.qualifiers.keys() and
                 (len(isoformMatch)==0 or
                  isoformMatch[0]=='1' or
                  isoformMatch[0]=='a')) or
                ('note' in f.qualifiers.keys() and
                 len(transcriptMatch)>0 and
                 (transcriptMatch[0]=='1' or
                  transcriptMatch[0]=='a' or
                  transcriptMatch[0]=='2'))):
                geneFound=True
                if (len(mRNA_exons)==0 and
                    len(several_mRNA_isoformsNums)==0 and
                    len(several_mRNA_isoformsWords)==0):
                    print('WARNING (7)! mRNA feature was not found in the GenBank file '+refDir+chrs[t]+'.gb for the gene '+t+'!')
                    mRNA_exons=f.location.parts
                    if args.includeNonCoding and t not in codons.keys():
                        writeRegions(resultFile,t,mRNA_exons,None,codons[geneNum],exons[geneNum],nucs[geneNum])
                        break
                elif (len(mRNA_exons)==0 and
                      (len(several_mRNA_isoformsNums)>0 or
                       len(several_mRNA_isoformsWords)>0)):
                    if len(several_mRNA_isoformsNums)>0:
                        mRNA_exons=sorted(several_mRNA_isoformsNums.items())[0][1]
                    else:
                        mRNA_exons=sorted(several_mRNA_isoformsWords.items())[0][1]
                    if args.includeNonCoding and t not in codons.keys():
                        writeRegions(resultFile,t,mRNA_exons,None,codons[geneNum],exons[geneNum],nucs[geneNum])
                        break
                writeRegions(resultFile,t,mRNA_exons,f.location.parts,codons[geneNum],exons[geneNum],nucs[geneNum])
                break
    if not geneFound:
        print('ERROR (8)! Gene with name "'+t+'" was not found in chromosome '+chrs[t])
        print(several_mRNA_isoformsWords)
        print(several_mRNA_isoformsNums)
        exit(8)
resultFile.close()
print('NGS-PrimerPlex finished!')
