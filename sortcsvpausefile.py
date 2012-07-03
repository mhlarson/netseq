import sys
import csv
import math

# Need to use the full pause file (i.e. with two lines for each pause)

def sortcsvpausefile(csvpausefilename):
  csvpausefile=open(csvpausefilename, 'rU')
  
  outfilename=csvpausefilename[:-4]

  #outfile_A=open(outfilename+"_A.csv",'wb')
  #outfile_C=open(outfilename+"_C.csv",'wb')
  #outfile_G=open(outfilename+"_G.csv",'wb')
  #outfile_T=open(outfilename+"_T.csv",'wb')
  
  A_counts=0
  C_counts=0
  G_counts=0
  T_counts=0
  
  csvpausestring=csvpausefile.readlines()
  for i in range(0,len(csvpausestring)):
    if csvpausestring[i][15]=='A':
      A_counts=A_counts+1
      #outfile_A.write(csvpausestring[i])
      #outfile_A.write(csvpausestring[i])  
    elif csvpausestring[i][15]=='C':
      C_counts=C_counts+1
      #outfile_C.write(csvpausestring[i])
      #outfile_C.write(csvpausestring[i]) 
    elif csvpausestring[i][15]=='G':
      G_counts=G_counts+1
      #outfile_G.write(csvpausestring[i])
      #outfile_G.write(csvpausestring[i]) 
    elif csvpausestring[i][15]=='T':
      T_counts=T_counts+1
      #outfile_T.write(csvpausestring[i])
      #outfile_T.write(csvpausestring[i]) 
      
    
  TotalCounts = (A_counts)+(C_counts)+(G_counts)+(T_counts)
  percentA=math.floor((float(A_counts)/TotalCounts)*100)
  percentC=math.floor((float(C_counts)/TotalCounts)*100)
  percentG=math.floor((float(G_counts)/TotalCounts)*100)
  percentT=math.floor((float(T_counts)/TotalCounts)*100)
    
  print 'percent A', percentA
  print 'percent C', percentC
  print 'percent G', percentG
  print 'percent T', percentT
  
  #outfile_A.close()
  #outfile_C.close()
  #outfile_G.close()
  #outfile_T.close()
  csvpausefile.close()
  
def main():
  sortcsvpausefile(sys.argv[1])

if __name__ == '__main__':
  main()