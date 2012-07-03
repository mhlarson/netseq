import sys
import csv

def cleancsvpausefile(csvpausefilename):
  csvpausefile=open(csvpausefilename, 'rU')
  outfile=open(csvpausefilename[:-4]+"_clean.csv",'wb')
  #outfile=open(csvpausefilename[:-4]+"_expression.csv",'wb')

  for line in csvpausefile:
    if line[0]!='>':
      outfile.write(line)
    
   #if line[0]=='>':
   #   row=line.strip('\n').strip('>').split('_')
   #   index=row[0]
   #   strand=row[1]
   #   pausefold=row[-1]
   #   outfile.write(index+'_'+strand+'_'+pausefold+'\t'+pausefold+'\n')

  outfile.close()
  csvpausefile.close()
  
def main():
  cleancsvpausefile(sys.argv[1])

if __name__ == '__main__':
  main()