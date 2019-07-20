import argparse,xlrd,os
import xlsxwriter as xls

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

parser=argparse.ArgumentParser(description='This script adds tags to primers for NGS')
parser.add_argument('--input','-in',dest='input',type=str,help='input XLS-file with designed primers')
##parser.add_argument('--with-additional-names','-names',dest='addNames',action='store_true',help='use this parameter, if input-file contains additional names that you gave your primers manually. They should be inserted before default names')
parser.add_argument('--tags-file','-tags',dest='tagsFile',type=str,help='text file with tags that we want to add to each primer. Default: "'+thisDir+'kplex_for_primers.txt"',default=thisDir+'kplex_for_primers.txt')
args=parser.parse_args()

try:
    file=open(args.tagsFile)
except FileNotFoundError:
    print('ERROR! File with tags was not found:',args.tagsFile)
    exit(0)
addSeq=[]
for string in file:
    if 'For ' in string:
        continue
    addSeq.append(string.replace('\n','').replace('\r',''))
file.close()

wb=xlrd.open_workbook(args.input)
ws=wb.sheet_by_index(0)
primerNames=[]
primers=[]
for i in range(ws.nrows):
    row=ws.row_values(i)
    if i==0 or row[0]=='':
        continue
    if row[1]!='':
        primerNames.append(row[3]+'_F')
        primers.append(addSeq[0]+row[1])
    if row[2]!='':
        primerNames.append(row[3]+'_R')
        primers.append(addSeq[1]+row[2])
# If there is also external primers
if len(wb.sheet_names())>1:
    ws=wb.sheet_by_index(1)
    for i in range(ws.nrows):
        row=ws.row_values(i)
        if i==0 or row[0]=='':
            continue
        if row[1]!='':
            primerNames.append(row[3]+'_F')
            primers.append(row[1])
        if row[2]!='':
            primerNames.append(row[3]+'_R')
            primers.append(row[2])
wbw=xls.Workbook(args.input[:-4]+'_with_seq.xls')
wsw=wbw.add_worksheet('Primers')
wsw.write_row(0,0,['Primer',"5' to 3' sequence"])
wsw.set_column(0,0,25)
wsw.set_column(1,1,60)
for i,(name,primer) in enumerate(zip(primerNames,primers)):
    wsw.write_row(i+1,0,[name,primer])
wbw.close()
