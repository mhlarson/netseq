import sys
import csv

def genelength(TSSfilename,outfilename):
  TSSfile=open(TSSfilename, 'rU')
  outfile=open(outfilename, 'wb')
  
  for line in TSSfile:
    start=int(line.split('\t')[1])
    stop=int(line.split('\t')[2])
    genelength=stop-start
    outfile.write(str(genelength)+'\n')
      
  TSSfile.close()
  outfile.close()
  
def main():
  genelength(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()