import sys
import subprocess
import random

def foldRNA(pausefilename,outfilename):
  pausefile=open(pausefilename, 'rU')
  outfile=open(outfilename,'wb')
   
  for line in pausefile:
    pausestring=line.strip('\n')
    pauselist=list(pausestring)
    random.shuffle(pauselist)
    pausestring=''.join(pauselist)
    currentpausefile=open('c:\UNAFold\currentpauseseq','w')
    currentpausefile.write(pausestring)
    currentpausefile.close()
    mfold=subprocess.call('perl c:\UNAFold\UNAFold.pl c:\UNAFold\currentpauseseq',stderr=subprocess.STDOUT,stdout=subprocess.PIPE,shell=True)
    deltaGfile=open('c:\UNAFold\currentpauseseq.dG','r')
    deltaGstring=deltaGfile.readlines()
    deltaGfile.close()
    deltaG=deltaGstring[1].split('\t')[1]
    outfile.write(pausestring+'\t'+deltaG+'\n')
    
  pausefile.close()
  outfile.close()
  
def main():
  foldRNA(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()