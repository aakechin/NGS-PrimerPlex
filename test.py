# This script performs some tests of NGS-PrimerPlex:
# - Extracting regions
# - Primer design with all possible features
#   + blast
#   + SNPs
#   + embedded
#   + draft primers
#   + draft primers with embedded
#   + embedded for internal primers
#   + anchored PCR
# - Adding adapters

import os
import argparse
import sys
import subprocess as sp

thisDir=os.path.dirname(os.path.realpath(__file__))+os.sep

# Section of reading arguments
parser=argparse.ArgumentParser(description='This script performs some tests of NGS-PrimerPlex')
parser.add_argument('--refDir','-ref',
                    dest='refDir',type=str,
                    help="directory with reference files. It is required only if you don't use it in docker",
                    required=False)
parser.add_argument('--reference-genome','-wgref',
                    dest='wholeGenomeRef',type=str,
                    help="file with INDEXED whole-genome reference sequence. "
                         "It is required only if you don't use it in docker",
                    required=False)
parser.add_argument('--dbsnp-vcf','-dbsnp',
                 dest='dbSnpVcfFile',type=str,
                 help='VCF-file (may be gzipped) with dbSNP variations. '
                      "It is required only if you don't use it in docker or "
                      "are not going to use this feature",
                 required=False)
args=parser.parse_args()

print('\n### Global testing of almost all functions of NGS-PrimerPlex started (11 tests)! ###')

print('\n# 1. Testing extraction of genome regions... #')
cmd=['python3',thisDir+'getGeneRegions.py',
     '-glf',thisDir+os.sep.join(['test',
                                 'test_gene.txt']),
     '-rf',thisDir+os.sep.join(['test',
                                'test_gene.regions.csv'])]
if args.refDir:
    cmd.extend(['-ref',args.refDir])
else:
    cmd.extend(['-ref','/hg19'])
if args.wholeGenomeRef:
    cmd.extend(['-wgref',args.wholeGenomeRef])
else:
    cmd.extend(['-wgref','/hg19/ucsc.hg19.fasta'])
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(1)
    
print('\n# 2. Testing primer design with blast, SNPs, and adapter sequences... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1']
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(2)

print('\n# 3. Testing converting to draft-primers... #')
cmd=['python3',thisDir+'convertToDraftFile.py',
     '-in',thisDir+os.sep.join(['test',
                                'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls']),
     '-out',thisDir+os.sep.join(['test',
                                 'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.draft.xls'])]
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
##    sys.stdout.write("\r"+str(c))
    out.append(str(c))
##    sys.stdout.flush()
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(3)

print('\n# 4. Testing primer design with draft-primers... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-draft',thisDir+os.sep.join(['test',
                                   'test_gene.regions_NGS_primerplex_all_draft_primers.xls'])]
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(4)

print('\n# 5. Testing primer design for embedded PCR... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-embedded']
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(5)

print('\n# 6. Testing primer design for embedded PCR with draft-primers... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-embedded',
     '-draft',thisDir+os.sep.join(['test',
                                   'test_gene.regions_NGS_primerplex_all_draft_primers.xls'])]
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(6)

print('\n# 7. Testing primer design for embedded PCR with given internal primers... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-embedded',
     '-primers',thisDir+os.sep.join(['test',
                                     'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls'])]
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(7)

print('\n# 8. Testing primer design for embedded anchored PCR... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions_anchored.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-embedded']
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(8)

##print('\n# 9. Testing extraction of genome regions for A. thaliana... #')
##cmd=['python3',thisDir+'getGeneRegions.py',
##     '-glf',thisDir+os.sep.join(['test',
##                                 'test_gene_arabidopsis.txt']),
##     '-rf',thisDir+os.sep.join(['test',
##                                'test_gene_arabidopsis.regions.csv'])]
##if args.refDir:
##    cmd.extend(['-ref',args.refDir])
##else:
##    cmd.extend(['-ref','/hg19'])
##if args.wholeGenomeRef:
##    cmd.extend(['-wgref',args.wholeGenomeRef])
##else:
##    cmd.extend(['-wgref','/hg19/ucsc.hg19.fasta'])
##process=sp.Popen(cmd,shell=False,
##                 stdout=sp.PIPE,stderr=sp.STDOUT,
##                 universal_newlines=True)
##out=[]
##for c in iter(process.stdout.readline,''):
##    out.append(str(c))
##if 'NGS-PrimerPlex finished!' in '\n'.join(out):
##    print('Finished successfully!')
##else:
##    print('ERROR! This test was finished with error!')
##    print('\n'.join(out))
##    exit(1)
##
##print('\n# 10. Testing primer design for A. thaliana... #')
##cmd=['python3',thisDir+'NGS_primerplex.py',
##     '-regions',thisDir+os.sep.join(['test',
##                                     'test_gene_arabidopsis.regions.csv']),
##     '-th','2','-run','test',
##     '-blast',
##     '-ad1','ctctctatgggcagtcggtgatt',
##     '-ad2','ctgcgtgtctccgactcag',
##     '-returnvariantsnum','1']
##if args.wholeGenomeRef:
##    cmd.extend(['-ref',args.wholeGenomeRef])
##else:
##    cmd.extend(['-ref','/hg19/ucsc.hg19.fasta'])
##if args.dbSnpVcfFile:
##    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
##elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
##    cmd.extend(['-snps','-dbsnp','/hg19/common_all_20180423_hg19.vcf.gz'])    
##process=sp.Popen(cmd,shell=False,
##                 stdout=sp.PIPE,stderr=sp.STDOUT,
##                 universal_newlines=True)
##out=[]
##for c in iter(process.stdout.readline,''):
##    out.append(str(c))
##if 'NGS-PrimerPlex finished!' in '\n'.join(out):
##    print('Finished successfully!')
##else:
##    print('ERROR! This test was finished with error!')
##    print('\n'.join(out))
##    exit(10)

print('\n# 9. Testing addition of adapter sequences... #')
cmd=['python3',thisDir+'addSeqToPrimers.py',
     '-in',thisDir+os.sep.join(['test',
                                'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls'])]
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(9)
