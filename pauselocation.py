import sys
import csv

def pauselocation(csvpausefilename,TSSfilename,outfilename,genecountfilename):
  csvpausefile=open(csvpausefilename, 'rU')
  TSSfile=open(TSSfilename, 'rU')
  outfile=open(outfilename, 'wb')
  genecountfile=open(genecountfilename, 'wb')
  
  genecount={}
  for i in range(0,10000):
    genecount[i]=0
  
  distfromstart_histo={}
  for i in range(0,10000):
    distfromstart_histo[i]=0
  
  gene_dict={}  
  for line in TSSfile:
    start=int(line.split('\t')[1])
    stop=int(line.split('\t')[2])
    genelength=stop-start
    gene=line.split('\t')[12]
    gene_dict[gene]=genelength
    for j in range(0,genelength):
      genecount[j]+=1
  
  for line in csvpausefile:
    if line[0]=='>':
      gene=line.strip('>').split('_')[0]
      distfromstart=int(line.split('_')[3])
      distfromstart_histo[distfromstart]+=float(1)
      
  for k in range (0,2000):
    distfromstart_histo[k]/=genecount[k]
    outfile.write(str(distfromstart_histo[k])+'\n')
    genecountfile.write(str(genecount[k])+'\n')
  
  TSSfile.close()
  csvpausefile.close()
  outfile.close()
  genecountfile.close()
  
def main():
  pauselocation(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

if __name__ == '__main__':
  main()