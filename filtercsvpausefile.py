import sys
import csv

def filtercsvpausefile(csvpausefilename,outfilename):
  csvpausefile=open(csvpausefilename, 'rU')
  outfile=open(outfilename,'wb')

  #csvpausestring=csvpausefile.readlines()
  #for i in range(0,len(csvpausestring),2):
  for line in csvpausefile:
    if line[40]=='G':
      outfile.write(line[0:35]+'\n')
      #outfile.write(csvpausestring[i+1])  
    
  outfile.close()
  csvpausefile.close()
  
def main():
  filtercsvpausefile(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()