import sys
import csv
import numpy
import math

from numpy import *

def PausePrediction(genekeyfilename,outfilename,wigstrand):
  CDSfile = open ('CDS_MG1655.mochiview.txt','rU')
  CDSstring = CDSfile.readlines()
  CDSlist=[]
    
  fnafile = open('c:\Python26\MG1655.fna', 'rU')
  
  fnastring=fnafile.read()
  fnastring=fnastring.replace('\n','')
  
  genekeyfile=open(genekeyfilename)
  genekeylist=[]
    
  for line in genekeyfile:
    fields=line.replace('\n','').split('\t')
    gene=fields[0]
    genekeylist.append(gene)
  
  if wigstrand=='-':
    fnastring=fnastring.replace('A','t').replace('C','g').replace('G','c').replace('T','a')
    fnastring=fnastring.replace('a','A').replace('c','C').replace('g','G').replace('t','T')
  
  outfile = open(outfilename,'wb')
  
  
  for j in range(0,len(CDSstring)):
    CDSrow = CDSstring[j].strip('\n').split('\t')
    gene = CDSrow[5]
    start = int(CDSrow[1])
    stop = int(CDSrow[2])
    strand = CDSrow[3]
    if gene in genekeylist:
      CDSlist.append((gene,start,stop,strand))
    
  countvar=0
  
  if wigstrand=='+':
    for line in CDSlist:
      if line[3]=='+':
        genename=line[0]
        start_index=line[1]
        stop_index=line[2]
        for index in range(start_index,stop_index+1):
          #if fnastring[index-2:index+2]=='TGCG' and fnastring[index-10:index-8]=='GG':
          if fnastring[index-2:index+2]=='TGCG' and (fnastring[index-9]=='G' or fnastring[index-10]=='G'):
            countvar+=1
            outfile.write(genename+'\t'+wigstrand+'\t'+str(index+1)+'\n')
            
  if wigstrand=='-':
    for line in CDSlist:
      if line[3]=='-':
        genename=line[0]
        start_index=line[1]
        stop_index=line[2]
        for index in range(start_index,stop_index+1):
          #if fnastring[index-2:index+2]=='GCGT' and fnastring[index+8:index+10]=='GG':
          if fnastring[index-2:index+2]=='GCGT' and (fnastring[index+8]=='G' or fnastring[index+9]=='G'):
            countvar+=1
            outfile.write(genename+'\t'+wigstrand+'\t'+str(index)+'\n')
            
  print countvar
  
  CDSfile.close()
  fnafile.close()
        
def main():
  PausePrediction(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()