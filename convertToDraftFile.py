# This script converts NGS-PrimerPlex XLS-file with primers
# into draft file

import argparse
import xlrd
import xlsxwriter as xls

parser=argparse.ArgumentParser(description='This script converts NGS-PrimerPlex '
                                           'XLS-file with primers into draft file')
parser.add_argument('--input','-in',
                    dest='inFile',type=str,
                    help='input XLS-file with primers designed by NGS-PrimerPlex')
parser.add_argument('--output','-out',
                    dest='outFile',type=str,
                    help='file for outputing draft file for NGS-PrimerPlex')
args=parser.parse_args()

# Stores names for worksheets
sheetNames=['Draft_Internal_Primers',
            'Draft_External_Primers']

wb=xlrd.open_workbook(args.inFile)
wbw=xls.Workbook(args.outFile)
for sheetNum in range(wb.nsheets):
    wsw=wbw.add_worksheet(sheetNames[sheetNum])
    wsw.write_row(0,0,['Primer_Pair','Left_Primer_Start','Left_Primer_Length',
                       'Right_Primer_End','Right_Primer_Length',
                       'Left_Primer_Tm','Right_Primer_Tm',
                       'Amplicon_Length','Primers_Score','Chrom'])
    ws=wb.sheet_by_index(sheetNum)
    for i in range(ws.nrows):
        row=ws.row_values(i)
        if i==0 or row[4]=='':
            continue
        newRow=['_'.join(row[1:3]),row[5],len(row[1]),
                row[6],len(row[2]),
                row[10],row[11],
                row[7],10,row[4]]
        wsw.write_row(i,0,newRow)
wbw.close()
print('NGS-PrimerPlex finished!')
