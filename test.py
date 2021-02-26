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
import xlrd
import subprocess as sp

thisDir=os.path.dirname(os.path.realpath(__file__))+os.sep

# Section of functions
def checkExcelFile(fileName,
                   expectedRowNums,
                   expectedColNums,
                   expectedSheetNums=1,
                   emptyCells=None,
                   badPrimers=None,
                   badFilter=None):
    # Checking the main output file
    wb=xlrd.open_workbook(fileName)
    for j in range(wb.nsheets):
        if j+1>expectedSheetNums:
            break
        ws=wb.sheet_by_index(j)
        if ((type(expectedRowNums[j])==list and
             ws.nrows not in list(range(expectedRowNums[j][0],
                                        expectedRowNums[j][-1]+1))) or
            (type(expectedRowNums[j])==int and
             ws.nrows!=expectedRowNums[j])):
            print('Finished successfully '
                  'but the output contains incorrect number of rows:')
            print(fileName)
            print('The index of the worksheet:',j+1)
            print('The expected number of rows:',expectedRowNums[j])
            print('The obtained number of rows:',ws.nrows)
            exit(11)
        for i in range(ws.nrows):
            row=ws.row_values(i)
            if len(row)<expectedColNums[j]:
                print('Finished successfully '
                      'but the output contains incorrect number of columns:')
                print(fileName)
                print('The index of the worksheet:',j+1)
                print('The expected number of columns:',expectedColNums[j])
                print('The obtained number of columns:',len(row))
                exit(12)
            if ('' in row and
                (emptyCells==None or
                 row.count('')>1 or
                 row.index('')!=emptyCells[i])):
                print('Finished successfully '
                      'but the output contains empty cells!')
                print(fileName)
                print('The index of the worksheet:',j+1)
                print('The index of row:',i+1)
                print('The index of cell:',row.index('')+1)
                print('The expected emppy cells (by rows):',emptyCells)
                print('Row:',row)
                exit(13)
            if (j==0 and
                badPrimers!=None and
                badPrimers==row[0]):
                print('Finished successfully '
                      'but the output contains primers that should have been filtered out!')
                print(fileName)
                print('The index of the worksheet:',j+1)
                print('The index of row:',i+1)
                print('Primers that should have been filtered out due to '+badFilter+':',badPrimers)
                exit(14)

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

print('\n### Global testing of almost all functions of NGS-PrimerPlex started (9 tests)! ###')

print('\n# 1. Testing extraction of genome regions... #')
cmd=['python3',thisDir+'getGeneRegions.py',
     '-glf',thisDir+os.sep.join(['test',
                                 'test_gene.txt']),
     '-rf',thisDir+os.sep.join(['test',
                                'test_gene.regions.csv'])]
if args.refDir:
    cmd.extend(['-ref',args.refDir])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19'])
if args.wholeGenomeRef:
    cmd.extend(['-wgref',args.wholeGenomeRef])
else:
    cmd.extend(['-wgref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking output
    file=open(thisDir+os.sep.join(['test',
                                   'test_gene.regions.csv']))
    expectedCols=[['KRAS','1','B','NotW'],
                  ['EGFR','1','B','NotW'],
                  ['EGFR','1','B','NotW']]
    i=0
    for string in file:
        cols=string.replace('\n','').split('\t')
        if cols[3:]!=expectedCols[i]:
            print('Finished successfully '
                  'but the output differs from the expected one!')
            exit(10)
        i+=1
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
     '-returnvariantsnum','1',
     '-primernum1','5',
     '-minampllen','145',
     '-optampllen','150',
     '-maxampllen','150']
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking the main output file
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls']),
                   [3],[18])
    # Checking file with draft primers
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_all_draft_primers.xls']),
                   [[30,40]],[10])
    if args.wholeGenomeRef:
        # Checking file with draft primers after blast
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_specificity.xls']),
                       [[30,40]],[10])
    if args.dbSnpVcfFile:
        # Checking file with draft primers after SNPs
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_SNPs.xls']),
                       [[30,40]],[10])
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
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking output
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.draft.xls']),
                   [3],[10])
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(3)

print('\n# 4. Testing primer design with draft-primers, checking for specificity and overlapping SNPs... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions.csv']),
     '-th','2','-run','test',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-primernum1','5',
     '-minampllen','145',
     '-optampllen','150',
     '-maxampllen','150',
     '-draft',thisDir+os.sep.join(['test',
                                   'draft_primers_for_checking_specificity_and_SNPs.xls'])]
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking the main output file
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls']),
                   [3],[18])
    # Checking file with draft primers
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_all_draft_primers.xls']),
                   [35],[10])
    if args.wholeGenomeRef:
        # Checking file with draft primers after blast
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_specificity.xls']),
                       [[34,35]],[10],
                       badPrimers='CCTCTCCCTCCCTCCAG_TGTGTTCCCGGACATAGTC',
                       badFilter='non-target hybridizations')
    if args.dbSnpVcfFile:
        # Checking file with draft primers after SNPs
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_SNPs.xls']),
                       [[33,35]],[10],
                       badPrimers='CTCACCTCCACCGTGCAG_CCCGTATCTCCCTTCCCTGATTA',
                       badFilter='overlapping SNPs')
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
     '-maxextampllen','170',
     '-primernum1','3',
     '-embedded']
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking the main output file
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls']),
                   [3,3],[18,18],2)
    # Checking file with draft primers
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_all_draft_primers.xls']),
                   [[9,13],[9,13]],[10,10],2)
    if args.wholeGenomeRef:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_specificity.xls']),
                       [[9,13],[9,13]],[10,10],2)
    if args.dbSnpVcfFile:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_SNPs.xls']),
                       [[9,13],[9,13]],[10,10],2)
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
     '-maxextampllen','170',
     '-primernum1','3',
     '-embedded',
     '-draft',thisDir+os.sep.join(['test',
                                   'test_gene.regions_NGS_primerplex_all_draft_primers.xls'])]
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking the main output file
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls']),
                   [3,3],[18,18],2)
    # Checking file with draft primers
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_all_draft_primers.xls']),
                   [[9,13],[10,19]],[10,10],2)
    if args.wholeGenomeRef:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_specificity.xls']),
                       [[9,13],[10,19]],[10,10],2)
    if args.dbSnpVcfFile:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_SNPs.xls']),
                       [[9,13],[10,19]],[10,10],2)
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
     '-maxextampllen','170',
     '-primernum1','3',
     '-embedded',
     '-primers',thisDir+os.sep.join(['test',
                                     'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls'])]
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking the main output file
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_test_primers_combination_1_info.xls']),
                   [3,3],[18,18],2)
    # Checking file with draft primers
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_NGS_primerplex_all_draft_primers.xls']),
                   [[9,13],[9,13]],[10,10],2)
    if args.wholeGenomeRef:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_specificity.xls']),
                       [[9,13],[9,13]],[10,10],2)
    if args.dbSnpVcfFile:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_NGS_primerplex_all_draft_primers_after_SNPs.xls']),
                       [[9,13],[9,13]],[10,10],2)
    print('Finished successfully!')
else:
    print('ERROR! This test was finished with error!')
    print('\n'.join(out))
    exit(7)

print('\n# 8. Testing primer design for embedded anchored PCR... #')
cmd=['python3',thisDir+'NGS_primerplex.py',
     '-regions',thisDir+os.sep.join(['test',
                                     'test_gene.regions_anchored.csv']),
     '-th','2','-run','test_anchored',
     '-blast','-snps',
     '-ad1','ctctctatgggcagtcggtgatt',
     '-ad2','ctgcgtgtctccgactcag',
     '-returnvariantsnum','1',
     '-maxextampllen','170',
     '-primernum1','3',
     '-embedded']
if args.wholeGenomeRef:
    cmd.extend(['-ref',args.wholeGenomeRef])
else:
    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
if args.dbSnpVcfFile:
    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
elif os.path.exists('/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'):
    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
process=sp.Popen(cmd,shell=False,
                 stdout=sp.PIPE,stderr=sp.STDOUT,
                 universal_newlines=True)
out=[]
for c in iter(process.stdout.readline,''):
    out.append(str(c))
if 'NGS-PrimerPlex finished!' in '\n'.join(out):
    # Checking the main output file
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_anchored_NGS_primerplex_test_anchored_primers_combination_1_info.xls']),
                   [3,3],[18,18],2,[0,1,2])
    # Checking file with draft primers
    checkExcelFile(thisDir+os.sep.join(['test',
                                        'test_gene.regions_anchored_NGS_primerplex_all_draft_primers.xls']),
                   [[10,18],[10,18]],[10,10],2)
    if args.wholeGenomeRef:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_anchored_NGS_primerplex_all_draft_primers_after_specificity.xls']),
                       [[10,18],[10,18]],[10,10],2)
    if args.dbSnpVcfFile:
        # Checking file with draft primers
        checkExcelFile(thisDir+os.sep.join(['test',
                                            'test_gene.regions_anchored_NGS_primerplex_all_draft_primers_after_SNPs.xls']),
                       [[10,18],[10,18]],[10,10],2)
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
##    cmd.extend(['-ref','/NGS-PrimerPlex/hg19'])
##if args.wholeGenomeRef:
##    cmd.extend(['-wgref',args.wholeGenomeRef])
##else:
##    cmd.extend(['-wgref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
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
##    cmd.extend(['-ref','/NGS-PrimerPlex/hg19/ucsc.hg19.fasta'])
##if args.dbSnpVcfFile:
##    cmd.extend(['-snps','-dbsnp',args.dbSnpVcfFile])
##elif os.path.exists('/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'):
##    cmd.extend(['-snps','-dbsnp','/NGS-PrimerPlex/hg19/common_all_20180423_hg19.vcf.gz'])    
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
