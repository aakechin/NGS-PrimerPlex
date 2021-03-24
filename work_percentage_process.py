# This script shows percent of performed work

# As input it takes:
# - number of done work (done)
# - number of whole work (wholeWork)

import sys

def showPercWork(done,wholeWork):
    percDoneWork=round((done/wholeWork)*100,2)
    sys.stdout.write("\r"+str(percDoneWork)+"%")
    sys.stdout.flush()
