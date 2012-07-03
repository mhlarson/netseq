import sys
import csv
import math

# 

def csvseqfilestats(csvpausefilename,outfilename):
  csvpausefile=open(csvpausefilename, 'rU')
  csvpausestring=csvpausefile.readlines()
  
  outfile = open(outfilename,'wb')
  
  percentA=0
  percentC=0
  percentG=0
  percentT=0
  
  percent_list=[]
  
  for i in range(0,len(csvpausestring)):
    A_counts=0
    C_counts=0
    G_counts=0
    T_counts=0
    for j in range(0,len(csvpausestring[i])):
      if csvpausestring[i][j]=='A':
        A_counts=A_counts+1
      if csvpausestring[i][j]=='C':
        C_counts=C_counts+1
      if csvpausestring[i][j]=='G':
        G_counts=G_counts+1
      if csvpausestring[i][j]=='T':
        T_counts=T_counts+1
    
    TotalCounts = (A_counts)+(C_counts)+(G_counts)+(T_counts)
    percentA=math.floor((float(A_counts)/TotalCounts)*100)
    percentC=math.floor((float(C_counts)/TotalCounts)*100)
    percentG=math.floor((float(G_counts)/TotalCounts)*100)
    percentT=math.floor((float(T_counts)/TotalCounts)*100)
    
    percent_list.append((i,percentA,percentC,percentG,percentT)) 
    
  csvlist = csv.writer(outfile, delimiter = '\t')
  for row in percent_list:
    csvlist.writerow(row)
    
  csvpausefile.close()
  outfile.close()
  
def main():
  csvseqfilestats(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()