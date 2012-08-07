import sys
import subprocess
import random

## Need to run the python command in the c:\UNAfold folder, otherwise program won't properly write folding energies into the output file ##

def foldRNA(pausefilename,outfilename):
  pausefile=open(pausefilename, 'rU') ## full version of longseq pause file ##
  outfile=open(outfilename,'wb')
  outfile.write('pause_fold'+'\t'+'upstream_sequence'+'\t'+'deltaG'+'\n')
  
  for line in pausefile:
    if line[0]=='>':
      pauserow=line.strip('>').strip('\n').split('_')
      pause_fold=pauserow[8]
      outfile.write(pause_fold+'\t')
    elif line[0]!='>':
      pausestring=line.strip('\n')[0:29]
      #pauselist=list(pausestring)
      #random.shuffle(pauselist)
      #pausestring=''.join(pauselist)
      currentpausefile=open('c:\UNAFold\currentpauseseq','w')
      currentpausefile.write(pausestring)
      currentpausefile.close()
      mfold=subprocess.call('perl c:\UNAFold\UNAFold.pl --constraints=c:\UNAFold\constraints.txt c:\UNAFold\currentpauseseq',stderr=subprocess.STDOUT,stdout=subprocess.PIPE,shell=True)
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