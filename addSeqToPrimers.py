import argparse,xlrd,os
import xlsxwriter as xls

thisDir=os.path.dirname(os.path.realpath(__file__))+'/'

parser=argparse.ArgumentParser(description='This script adds tags to primers for NGS')
parser.add_argument('--input','-in',
                    dest='input',type=str,
                    help='input XLS-file with designed primers')
parser.add_argument('--tags-file','-tags',
                    dest='tagsFile',type=str,
                    help='text file with tags that we want to add to each primer. Default: "'+thisDir+'kplex_for_primers.txt"',
                    default=thisDir+'kplex_for_primers.txt')
parser.add_argument('--add-cells','-cells',
                    dest='addCells',action='store_true',
                    help='use this argument, if you want to add cell ID for the synthesis in plate')
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
wbw=xls.Workbook(args.input[:-4]+'_with_seq.xls')
wswAll=wbw.add_worksheet('All_Primers')
wsws={'All_Primers':wswAll}
wswRowNum={'All_Primers':1}
wswLetters={'All_Primers':ord('A')}
wswNums={'All_Primers':1}
wswAll.write_row(0,0,['Primer',"5' to 3' sequence"])
wswAll.set_column(0,0,25)
wswAll.set_column(1,1,60)
for i in range(ws.nrows):
    row=ws.row_values(i)
    if i==0 or row[0]=='':
        continue
    if row[17]!='':
        multNum=str(int(row[17]))
        if multNum not in wsws.keys():
            # Check if all previous multiplexes are absent, too
            for i in range(1,int(multNum)):
                if str(i) not in wsws.keys():
                    wsw=wbw.add_worksheet(str(i))
                    wsws[str(i)]=wsw
                    wswRowNum[str(i)]=1
                    wswLetters[str(i)]=ord('A')
                    wswNums[str(i)]=1
                    wsw.write_row(0,0,['Primer',"5' to 3' sequence"])
                    wsw.set_column(0,0,25)
                    wsw.set_column(1,1,60)
            wsw=wbw.add_worksheet(multNum)
            wsws[multNum]=wsw
            wswRowNum[multNum]=1
            wswLetters[multNum]=ord('A')
            wswNums[multNum]=1
            wsw.write_row(0,0,['Primer',"5' to 3' sequence"])
            wsw.set_column(0,0,25)
            wsw.set_column(1,1,60)
        else:
            wsw=wsws[multNum]
    if row[1]!='':
        wswAll.write_row(wswRowNum['All_Primers'],0,[row[3]+'_F',
                                                     addSeq[0]+row[1]])
        if args.addCells:
            wswAll.write(wswRowNum['All_Primers'],2,chr(wswLetters['All_Primers'])+str(wswNums['All_Primers']))
        wswRowNum['All_Primers']+=1
        if row[17]!='':
            wsw.write_row(wswRowNum[multNum],0,[row[3]+'_F',
                                                     addSeq[0]+row[1]])
            if args.addCells:
                wsw.write(wswRowNum[multNum],2,chr(wswLetters[multNum])+str(wswNums[multNum]))
            wswRowNum[multNum]+=1
    if row[2]!='':
        wswAll.write_row(wswRowNum['All_Primers'],0,[row[3]+'_R',
                                                     addSeq[1]+row[2]])
        if args.addCells:
            wswAll.write(wswRowNum['All_Primers'],2,chr(wswLetters['All_Primers'])+str(wswNums['All_Primers']+1))
        wswRowNum['All_Primers']+=1
        if row[17]!='':
            wsw.write_row(wswRowNum[multNum],0,[row[3]+'_R',
                                                     addSeq[1]+row[2]])
            if args.addCells:
                wsw.write(wswRowNum[multNum],2,chr(wswLetters[multNum])+str(wswNums[multNum]+1))
            wswRowNum[multNum]+=1
    if args.addCells:
        if wswLetters['All_Primers']==ord('H'):
            wswLetters['All_Primers']=ord('A')
            if wswNums['All_Primers']==11:
                wswNums['All_Primers']=1
            else:
                wswNums['All_Primers']+=2
        else:
            wswLetters['All_Primers']+=1
        if row[17]!='':
            if wswLetters[multNum]==ord('H'):
                wswLetters[multNum]=ord('A')
                if wswNums[multNum]==11:
                    wswNums[multNum]=1
                else:
                    wswNums[multNum]+=2
            else:
                wswLetters[multNum]+=1
# If there is also external primers
if len(wb.sheet_names())>1:
    ws=wb.sheet_by_index(1)
    for i in range(ws.nrows):
        row=ws.row_values(i)
        if i==0 or row[0]=='':
            continue
        if row[18]!='':
            multNum=str(int(row[18]))
        if row[1]!='':
            wswAll.write_row(wswRowNum['All_Primers'],0,[row[3]+'_F',
                                                         addSeq[0]+row[1]])
            if args.addCells:
                wswAll.write(wswRowNum['All_Primers'],2,chr(wswLetters['All_Primers'])+str(wswNums['All_Primers']))
            wswRowNum['All_Primers']+=1
            if row[18]!='':
                wsw.write_row(wswRowNum[multNum],0,[row[3]+'_F',
                                                         addSeq[0]+row[1]])
                if args.addCells:
                    wsw.write(wswRowNum[multNum],2,chr(wswLetters[multNum])+str(wswNums[multNum]))
                wswRowNum[multNum]+=1
        if row[2]!='':
            wswAll.write_row(wswRowNum['All_Primers'],0,[row[3]+'_R',
                                                         addSeq[1]+row[2]])
            if args.addCells:
                wswAll.write(wswRowNum['All_Primers'],2,chr(wswLetters['All_Primers'])+str(wswNums['All_Primers']+1))
            wswRowNum['All_Primers']+=1
            if row[18]!='':
                wsw.write_row(wswRowNum[multNum],0,[row[3]+'_R',
                                                         addSeq[1]+row[2]])
                if args.addCells:
                    wsw.write(wswRowNum[multNum],2,chr(wswLetters[multNum])+str(wswNums[multNum]+1))
                wswRowNum[multNum]+=1
        if args.addCells:
            if wswLetters['All_Primers']==ord('H'):
                wswLetters['All_Primers']=ord('A')
                if wswNums['All_Primers']==11:
                    wswNums['All_Primers']=1
                else:
                    wswNums['All_Primers']+=2
            else:
                wswLetters['All_Primers']+=1
            if row[18]!='':
                if wswLetters[multNum]==ord('H'):
                    wswLetters[multNum]=ord('A')
                    if wswNums[multNum]==11:
                        wswNums[multNum]=1
                    else:
                        wswNums[multNum]+=2
                else:
                    wswLetters[multNum]+=1

wbw.close()
print('NGS-PrimerPlex finished!')
