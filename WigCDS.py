import sys
import csv
import numpy
import math

# Runs a cross-correlation on two wig files in which each index is represented, even if there are zero counts at a given index.
# Riboseqfilename is moved forward/downstream across Netseqfilename over a given range. In bacteria, a transcriptional pause
# at a given sequence should lead to a ribosomal pause upstream of the transcriptional pause site. This program determines the 
# downstream shift in riboseqfilename that gives the maximal correlation. For plus strand, use as directed. For minus strand, reverse
# order of the input files. 

def wigCDS(wigfilename,outwigfilename,wigstrand):
  CDSfile = open ('CDS_MG1655.mochiview.txt','rU')
  CDSstring = CDSfile.readlines()
  CDSdict={}
  CDSlist=[]
  
  for j in range(0,len(CDSstring)):
    CDSrow = CDSstring[j].strip('\n').split('\t')
    gene = CDSrow[5]
    start = int(CDSrow[1])
    stop = int(CDSrow[2])
    strand = CDSrow[3]
    CDSlist.append((gene,start,stop,strand))
  
  if wigstrand=='plus':
    for line in CDSlist:
      if line[3]=='+':
        for index in range(line[1],line[2]+1):
          CDSdict[index]=line[0]
          
  elif wigstrand=='minus':
    for line in CDSlist:
      if line[3]=='-':
        for index in range(line[1],line[2]+1):
          CDSdict[index]=line[0]
      
  wigfile = open(wigfilename, 'rU')
  wigstring = wigfile.readlines()
  wiglist = []
  outwigfile = open(outwigfilename, 'wb')
  
  for j in range(len(wigstring)-1):
    wigrow = wigstring[j].strip('\n').split(' ')
    index = int(wigrow[0])
    count = float(wigrow[1])
    if index in CDSdict:
      wiglist.append((index,count))
          
  csvlist = csv.writer(outwigfile, delimiter = ' ')
  for row in wiglist:
    csvlist.writerow(row)
  
  wigfile.close()
  outwigfile.close()
      
def main():
  wigCDS(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()