import sys
import csv
import numpy
import math

def ExtractPauseFold(pause_filename,outfilename):
  pause_file = open(pause_filename, 'rU')
  pause_list=[]
  pause_fold_list=[]
  outfile = open(outfilename,'wb')
  
  j=0
  for line in pause_file:
    if line[0]=='>':
      pause_row = line.strip('\n').split('_')
      pause_index = int(pause_row[4])
      pause_fold = float(pause_row[-1])
      pause_list.append((pause_index,pause_fold))
      j+=1
  
  for row in pause_list:
    pause_fold_list.append((row[0],row[1]))
  
  csvlist = csv.writer(outfile, delimiter = '\t')
  for row in pause_fold_list:
    csvlist.writerow(row)
  
  pause_file.close()
  outfile.close()

def main():
  ExtractPauseFold(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()