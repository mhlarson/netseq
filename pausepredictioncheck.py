import sys
import csv

def pausepredictioncheck(consensussitesfilename,pausefilename,outfilename):
  
  consensussitesfile=open(consensussitesfilename, 'rU') #This is consensus pause sites for only well-transcribed genes.
  pausefile=open(pausefilename, 'rU')                   #Pause file with two lines for each pause
  outfile = open(outfilename,'wb')
  
  consensus_dict={}
  
  for line in consensussitesfile:
    row=line.strip('\n').split('\t')
    consensus_index=row[2]
    consensus_strand=row[1]
    consensus_gene=row[0]
    consensus_dict[consensus_index]=consensus_gene
  
  for line in pausefile:
    if line[0]=='>':
      row=line.strip('\n').strip('>').split('_')
      index=row[0]
      strand=row[1]
      pausefold=row[-1]
      if index in consensus_dict:
        outfile.write(index+'_'+strand+'_'+pausefold+'\t'+pausefold+'\n')

  consensussitesfile.close()
  outfile.close()
  pausefile.close()
  
def main():
  pausepredictioncheck(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()