import sys
import csv
import numpy
import math

def PauseFilter(pauseA_filename,pauseB_filename,outpausefilename,strand):
  pauseA_file = open(pauseA_filename, 'rU')
  pauseA_dict = {}
  
  pauseB_file = open(pauseB_filename, 'rU')
  pauseB_dict = {}
  
  pauseAB_list=[]
  pauseAnoB_list=[]
  pauseBnoA_list=[]
  backtrackpause_list=[]
   
  outpausefileAB = open(outpausefilename+'_greA_and_mRNA.csv','wb')
  outpausefileAnoB = open(outpausefilename+'_greA_nomRNA.csv','wb')
  outpausefileBnoA = open(outpausefilename+'_mRNA_nogreA.csv','wb')
  #outpausefileBacktrack=open(outpausefilename+'_bactrackpauses.csv','wb')
  
  fnafile = open('c:\Python26\MG1655.fna', 'rU')
  
  fnastring=fnafile.read()
  fnastring=fnastring.replace('\n','')
  
  if strand=='-':
    fnastring=fnastring.replace('A','t').replace('C','g').replace('G','c').replace('T','a')
    fnastring=fnastring.replace('a','A').replace('c','C').replace('g','G').replace('t','T')
    
  for line in pauseA_file:
    if line[0]=='>':
      pauseA_row = line.strip('\n').split('_')
      pauseA_index = int(pauseA_row[4])
      pauseA_fold = float(pauseA_row[-1])
      pauseA_strand=pauseA_row[1]
      if pauseA_strand==strand:
        pauseA_dict[pauseA_index]=pauseA_fold
                
  for line in pauseB_file:
    if line[0]=='>':
      pauseB_row = line.strip('\n').split('_')
      pauseB_index = int(pauseB_row[4])
      pauseB_fold = float(pauseB_row[-1])
      pauseB_strand=pauseB_row[1]
      if pauseB_strand==strand:
        pauseB_dict[pauseB_index]=pauseB_fold
      
  for row in pauseA_dict:
    if row in pauseB_dict:
      outpausefileAB.write('>'+str(row)+'_'+strand+'_'+str(pauseA_dict[row])+'\n')
      outpausefileAB.write(fnastring[row-15:row+6]+'\n')
    elif row not in pauseB_dict and row+1 not in pauseB_dict and row+2 not in pauseB_dict and row+3 not in pauseB_dict:
      outpausefileAnoB.write('>'+str(row)+'_'+strand+'_'+str(pauseA_dict[row])+'\n')
      outpausefileAnoB.write(fnastring[row-15:row+6]+'\n')
  
  #for row in pauseA_dict:
  #  if row not in pauseB_dict and row+2 in pauseB_dict or row+3 in pauseB_dict:
  #    backtrackpause_list.append((row,fnastring[row-15:row+6]))
    
  for row in pauseB_dict:
    if row not in pauseA_dict:
      outpausefileBnoA.write('>'+str(row)+'_'+strand+'_'+str(pauseB_dict[row])+'\n')
      outpausefileBnoA.write(fnastring[row-15:row+6]+'\n')
      
  #csvlist = csv.writer(outpausefileAB, delimiter = '\t')
  #for row in pauseAB_list:
  #  csvlist.writerow(row)
    
  #csvlist = csv.writer(outpausefileAnoB, delimiter = '\t')
  #for row in pauseAnoB_list:
  #  csvlist.writerow(row)
  
  #csvlist = csv.writer(outpausefileBnoA, delimiter = '\t')
  #for row in pauseBnoA_list:
  #  csvlist.writerow(row)
    
  #csvlist = csv.writer(outpausefileBacktrack, delimiter = '\t')
  #for row in backtrackpause_list:
  #  csvlist.writerow(row)
  
  pauseA_file.close()
  pauseB_file.close()
  outpausefileAB.close()
  outpausefileAnoB.close()
  outpausefileBnoA.close()
  #outpausefileBacktrack.close()

def main():
  PauseFilter(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

if __name__ == '__main__':
  main()