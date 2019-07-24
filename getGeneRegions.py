# This script takes names of genes and numbers of their exons or positions in CDS
# and makes region file for NGS-PrimerPlex program

# Section of importing modules
from Bio import SeqIO
from Bio import Entrez as en
import copy
import argparse
import os
import re
import gzip

# Global variables
global thisDir
thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

# Section of functions
def getChrNum(gName,org):
    h=en.esearch(term=gName+'[Gene Name]',db='gene')
    r=en.read(h)
    for i in r['IdList']:
        h=en.esummary(db='gene',id=i)
        r2=en.read(h)
        fName=r2['DocumentSummarySet']['DocumentSummary'][0]['NomenclatureSymbol']
        fOrg=r2['DocumentSummarySet']['DocumentSummary'][0]['Organism']
        fSyn=r2['DocumentSummarySet']['DocumentSummary'][0]['NomenclatureSymbol'].split(', ')
        fChr=r2['DocumentSummarySet']['DocumentSummary'][0]['Chromosome']
        if fOrg['CommonName']==org and (fName==gName or fName in fSyn):
            return(fChr)

def writeRegions(resultFile,targetGene,pro_mRNA,cds=None,codons={},exons={},nucs={}):
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
    if targetGene not in codons.keys() and targetGene not in exons.keys():
        for exonNum,subf in mRNA.items():
            if exonNum+1!=1 and exonNum+1!=len(mRNA):
                resultFile.write('chr'+chrs[targetGene]+'\t'+str(subf[0]+1-intronSize)+'\t'+str(subf[1]+intronSize)+'\t'+targetGene+'\n')
            elif (exonNum+1==1 and subf[2]>0) or (exonNum+1==len(mRNA) and subf[2]<0):
                resultFile.write('chr'+chrs[targetGene]+'\t'+str(subf[0]+1)+'\t'+str(subf[1]+intronSize)+'\t'+targetGene+'\n')
            elif (exonNum+1==1 and subf[2]<0) or (exonNum+1==len(mRNA) and subf[2]>0):
                resultFile.write('chr'+chrs[targetGene]+'\t'+str(subf[0]+1-intronSize)+'\t'+str(subf[1])+'\t'+targetGene+'\n')
            else:
                print('ERROR: Unknown variant of exon number and gene strand!')
                print(mRNA)
                print(exonNum+1,subf[2])
                exit(1)
    elif targetGene in codons.keys():
        exLen=0
        exLens={}
        # If the 1st coding exon is no the 1st exon
        for exonNum in range(min(mRNA.keys())):
            exLens[exonNum]=0
        for exonNum,subf in sorted(mRNA.items()):
            exLen+=subf[1]-subf[0]
            exLens[exonNum]=exLen
        for reg in nucs[targetGene]:
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
                print('ERROR: exon for region '+str(reg)+' of gene '+targetGene+' was not determined correctly')
                exit(0)
##                        print(exNum1,exNum2,fromStart1,fromStart2)
##                        print(f.location.parts[exNum1].nofuzzy_end,f.location.parts[exNum2].nofuzzy_start)
            # If exons for the start and for the end of region are equal
            if exNum1==exNum2:
                if list(mRNA.values())[0][2]==-1:
                    try:
                        resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum2][1]-fromStart2)+'\t'+str(mRNA[exNum1][1]-fromStart1)+'\t'+targetGene+'\n')
                    except KeyError:
                        print('ERROR!',exNum1,exNum2)
                        print(mRNA)
                        exit(1)
                else:
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum1][0]+1+fromStart1)+'\t'+str(mRNA[exNum2][0]+1+fromStart2)+'\t'+targetGene+'\n')
            elif exNum2-exNum1==1:
                if list(mRNA.values())[0][2]==-1:
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum1][0]+1)+'\t'+str(mRNA[exNum1][1]-fromStart1)+'\t'+targetGene+'\n')
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum2][1]-fromStart2)+'\t'+str(mRNA[exNum2][1])+'\t'+targetGene+'\n')
                else:
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum1][0]+1+fromStart1)+'\t'+str(mRNA[exNum1][1])+'\t'+targetGene+'\n')
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum2][0]+1)+'\t'+str(mRNA[exNum2][0]+1+fromStart2)+'\t'+targetGene+'\n')
            else:
                if list(mRNA.values())[0][2]==-1:
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum1][0]+1)+'\t'+str(mRNA[exNum1][1]-fromStart1)+'\t'+targetGene+'\n')
                    for j in range(exNum1+1,exNum2):
                        resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[j][0]+1)+'\t'+str(mRNA[j][1])+'\t'+targetGene+'\n')
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum2][1]-fromStart2)+'\t'+str(mRNA[exNum2][1])+'\t'+targetGene+'\n')
                else:
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum1][0]+1+fromStart1)+'\t'+str(mRNA[exNum1][1])+'\t'+targetGene+'\n')
                    for j in range(exNum1+1,exNum2):
                        resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[j][0]+1)+'\t'+str(mRNA[j][1])+'\t'+targetGene+'\n')
                    resultFile.write('chr'+chrs[targetGene]+'\t'+str(mRNA[exNum2][0]+1)+'\t'+str(mRNA[exNum2][0]+1+fromStart2)+'\t'+targetGene+'\n')
    elif targetGene in exons.keys():
        for exRange in exons[targetGene]:
##            print(mRNA)
            for exonNum in range(exRange[0]-1,exRange[1]):
                if exonNum not in mRNA.keys():
                    print('ERROR! Exon '+str(exonNum)+' includes noncoding sequences or does not exist.')
                    print('In the first case, use parameter -noncoding. In the second one, correct the mistake')
                    exit(1)
                subf=mRNA[exonNum]
                resultFile.write('chr'+chrs[targetGene]+'\t'+str(subf[0]+1)+'\t'+str(subf[1])+'\t'+'_'.join([targetGene,str(exonNum+1)])+'\n')

# Section of reading arguments
parser=argparse.ArgumentParser(description='This script takes names of genes '
                               'and numbers of their exons or positions in CDS '
                               'and makes region file for hi-plex program')
parser.add_argument('--geneListFile','-glf',type=str,
                    help='file with list of genes. Format is: GENE EXONS CODONS',required=True)
parser.add_argument('--refDir','-ref',type=str,
                    help='directory with reference files',required=True)
parser.add_argument('--organism','-org',type=str,
                    help='common name of organism (human, sheep, etc). Default: human',default='human')
parser.add_argument('--resultFile','-rf',type=str,
                    help='file for results',required=True)
parser.add_argument('--intron-nucleotides','-intron',type=int,dest='intronSize',
                    help='number of nucleotides from intron to take. Default: 2',default=2)
parser.add_argument('--include-noncoding','-noncoding',dest='includeNonCoding',action='store_true',
                    help="use this parameter, if you want to include 5'- and 3'-non-coding regions of mRNA")
parser.add_argument('--email','-email',type=str,
                    help='e-mail for getting chromosome number '
                         'from the ENTREZ database for the user-defined genes.'
                         ' Default: info@gmail.com',
                    required=False,default='info@gmail.com')
args=parser.parse_args()

en.email=args.email
refDir=args.refDir
org=args.organism
intronSize=args.intronSize

# Read input file
print('Reading input file and getting chromosome numbers...')
targetGenes=[]
chrs={}
codons={}
exons={}
geneListFile=open(args.geneListFile)
for string in geneListFile:
    if 'Gene\tExons' in string:
        continue
    cols=string[:-1].split('\t')
    print(cols[0])
    targetGenes.append(cols[0])
    chrs[cols[0]]=getChrNum(cols[0],org)
    if chrs[cols[0]]==None:
        print('ERROR: Gene '+cols[0]+' was not found in the Gene database of NCBI!')
        exit(1)
    if not (len(cols)<2 or cols[1]==''):
        ex=cols[1].split(',')
        exons[cols[0]]=[]
        for e in ex:
            if '-' in e:
                ep=e.split('-')
                exons[cols[0]].append([int(ep[0]),int(ep[1])])
            else:
                exons[cols[0]].append([int(e),int(e)])
    if not (len(cols)<3 or cols[2]==''):
        ex=cols[2].split(',')
        codons[cols[0]]=[]
        for e in ex:
            if '-' in e:
                ep=e.split('-')
                codons[cols[0]].append([int(ep[0]),int(ep[1])])
            else:
                codons[cols[0]].append([int(e),int(e)])
geneListFile.close()
print('Done')

# Library nucs contains ranges of nucleotides for which we need to line up primers
# we get it from codons list
nucs=copy.deepcopy(codons)
for key,value in nucs.items():
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
for t in targetGenes:
    print(t)
    try:
        data=SeqIO.read(refDir+'chr'+chrs[t]+'.gb','genbank')
    except FileNotFoundError:
        try:
            data=SeqIO.read(gzip.open(refDir+'chr'+chrs[t]+'.gb.gz','rt'),'genbank')
        except:
            print('ERROR: GenBank file for chromosome '+chrs[t]+' is absent in the reference directory '+refDir+'!')
            exit(1)
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
##                print('mRNA1')
                mRNA_exons=f.location.parts
                if args.includeNonCoding and t not in codons.keys():
                    geneFound=True
                    writeRegions(resultFile,t,mRNA_exons,None,codons,exons,nucs)
                    break
            elif (len(transcriptPat.findall(f.qualifiers['product'][0]))==0 or
                transcriptPat.findall(f.qualifiers['product'][0])[0]=='1' or
                transcriptPat.findall(f.qualifiers['product'][0])[0]=='a'):
                mRNA_exons=f.location.parts
##                print('mRNA2')
                if args.includeNonCoding and t not in codons.keys():
                    geneFound=True
                    writeRegions(resultFile,t,mRNA_exons,None,codons,exons,nucs)
                    break
            elif len(transcriptPat.findall(f.qualifiers['product'][0]))>0:
##                print('mRNA3')
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
##            print('CDS0',isoformMatch,f.qualifiers['product'][0],transcriptMatch,f.qualifiers['note'][0])
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
                  transcriptMatch[0]=='a'))):
                geneFound=True
##                print('CDS1')
                if (len(mRNA_exons)==0 and
                    len(several_mRNA_isoformsNums)==0 and
                    len(several_mRNA_isoformsWords)==0):
                    print('ERROR: mRNA feature was not found in the GenBank file '+refDir+'chr'+chrs[t]+'.gb for the gene '+t+'!')
                    exit(1)
                elif (len(mRNA_exons)==0 and
                      (len(several_mRNA_isoformsNums)>0 or
                       len(several_mRNA_isoformsWords)>0)):
                    if len(several_mRNA_isoformsNums)>0:
                        mRNA_exons=sorted(several_mRNA_isoformsNums.items())[0][1]
                    else:
                        mRNA_exons=sorted(several_mRNA_isoformsWords.items())[0][1]
                    if args.includeNonCoding and t not in codons.keys():
                        writeRegions(resultFile,t,mRNA_exons,None,codons,exons,nucs)
                        break
                writeRegions(resultFile,t,mRNA_exons,f.location.parts,codons,exons,nucs)
                break
    if not geneFound:
        print('ERROR: gene with name "'+t+'" was not found in chromosome '+chrs[t])
        exit(0)
resultFile.close()
print('Done')
