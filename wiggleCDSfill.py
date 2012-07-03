import sys
import csv
import numpy
import math

def wiggleCDSfill(wigfilename,outwigfilename):
  wigfile = open(wigfilename, 'rU')
  wigstring = wigfile.readlines()
  wigdict = {}
  outwigfile = open(outwigfilename+".wig", 'wb')
  
  CDSfile = open ('CDS_MG1655.mochiview.txt','rU')
  CDSstring = CDSfile.readlines()
  CDSlist=[]
  
  geneReads = []
  
  countlis=[]
  
  for j in range(0,len(CDSstring)):
    CDSrow = CDSstring[j].strip('\n').split('\t')
    gene = CDSrow[5]
    start = int(CDSrow[1])
    stop = int(CDSrow[2])
    strand = CDSrow[3]
    CDSlist.append((gene,start,stop,strand))
  
  wigstrand='plus'
  
  if wigstrand=='plus':
    for line in CDSlist:
      if line[3]=='+':
        genename=line[0]
        start_index=line[1]
        stop_index=line[2]
        for index in range(start_index,stop_index+1):
          wigdict[index]=0
  
  for j in range(len(wigstring)-1):
    wigrow = wigstring[j].strip('\n').split(' ')
    index = int(wigrow[0])
    count = float(wigrow[1])
    if index in wigdict:
      wigdict[index]=count

  wiglist=sorted(wigdict.items())
    
  
  
  csvlist = csv.writer(outwigfile , delimiter = '\t')
  for row in wiglist:
    csvlist.writerow(row)
  
  wigfile.close()
  outwigfile.close()
  CDSfile.close()

def main():
  wiggleCDSfill(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()