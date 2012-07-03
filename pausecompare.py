import sys
import csv
import numpy
import math

def PauseCompare(genekeyfilename,pauseA_filename,pauseB_filename,outpausefilename,strand):
  
  genekeyfile=open(genekeyfilename)
  genekeydict={}
  
  for line in genekeyfile:
    fields=line.replace('\n','').split('\t')
    gene=fields[0]
    transcriptionlevel=fields[2]
    genekeydict[gene]=transcriptionlevel
  
  
  pauseA_file = open(pauseA_filename, 'rU')
  pauseA_dict = {}
  
  pauseB_file = open(pauseB_filename, 'rU')
  pauseB_list = []
  
  pauseBA_list=[]
  pauseBnoA_list=[]
  
  outpausefileBA = open(outpausefilename+'BA.csv','wb')
  outpausefileBnoA = open(outpausefilename+'BnoA.csv','wb')
  
  fnafile = open('c:\Python26\MG1655.fna', 'rU')
  
  fnastring=fnafile.read()
  fnastring=fnastring.replace('\n','')
  
  if strand=='-':
    fnastring=fnastring.replace('A','t').replace('C','g').replace('G','c').replace('T','a')
    fnastring=fnastring.replace('a','A').replace('c','C').replace('g','G').replace('t','T')
    
  for line in pauseA_file:
    if line[0]=='>':
      pauseA_row = line.strip('\n').strip('>').split('_')
      pauseA_gene=pauseA_row[0]
      pauseA_index = int(pauseA_row[4])
      pauseA_fold = float(pauseA_row[-1])
      pauseA_strand=pauseA_row[1]
      if pauseA_strand==strand:
        pauseA_dict[pauseA_index]=pauseA_fold
                
  for line in pauseB_file:
    if line[0]=='>':
      pauseB_row = line.strip('\n').strip('>').split('_')
      pauseB_gene=pauseB_row[0]
      pauseB_index = int(pauseB_row[4])
      pauseB_fold = float(pauseB_row[-1])
      pauseB_strand=pauseB_row[1]
      if pauseB_strand==strand:
        pauseB_list.append((pauseB_index,pauseB_fold,pauseB_gene))
      
  for row in pauseB_list:
    if row[2] in genekeydict and row[0] in pauseA_dict and pauseA_dict[row[0]]!=0:
      pauseBA_list.append((row[0],row[1]/pauseA_dict[row[0]],fnastring[row[0]-15:row[0]+6]))
    elif row[2] in genekeydict:
      pauseBnoA_list.append((row[0],fnastring[row[0]-15:row[0]+6]))
      
  csvlist = csv.writer(outpausefileBA, delimiter = '\t')
  for row in pauseBA_list:
    csvlist.writerow(row)
    
  csvlist = csv.writer(outpausefileBnoA, delimiter = '\t')
  for row in pauseBnoA_list:
    csvlist.writerow(row)
  
  pauseA_file.close()
  pauseB_file.close()
  outpausefileBA.close()
  outpausefileBnoA.close()

def main():
  PauseCompare(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

if __name__ == '__main__':
  main()