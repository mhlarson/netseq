import sys
import csv
import numpy
import math

from numpy import *

## Looks for consensus pause sequences in well-expressed genes. From that list, determines whether a pause occurs there or not. ##

def ConsensusPauseFilter(pausefilename,genekeyfilename,out_nopausefilename,out_pausefilename):
  CDSfile = open ('CDS_MG1655.mochiview.txt','rU')
  CDSstring = CDSfile.readlines()
  CDSlist=[]
    
  fnafile = open('c:\Python26\MG1655.fna', 'rU')
  
  fnastring=fnafile.read()
  fnastring=fnastring.replace('\n','')
  
  pausefile=open(pausefilename, 'rU')  ## Full pause file ##
  pause_dict={}
  
  genekeyfile=open(genekeyfilename)
  genekeydict={}
  
  for line in genekeyfile:
    fields=line.replace('\n','').split('\t')
    gene=fields[0]
    transcriptionlevel=fields[2]
    genekeydict[gene]=transcriptionlevel
  
  j=0
  for line in pausefile:
    if line[0]=='>':
      pause_row = line.strip('\n').split('_')
      pause_index = int(pause_row[4])
      pause_fold = float(pause_row[-1])
      pause_dict[pause_index]=pause_fold
      j+=1
  
  #if wigstrand=='-':
  #  fnastring=fnastring.replace('A','t').replace('C','g').replace('G','c').replace('T','a')
  #  fnastring=fnastring.replace('a','A').replace('c','C').replace('g','G').replace('t','T')
  
  out_nopausefile = open(out_nopausefilename,'wb')
  out_pausefile = open(out_pausefilename,'wb')
  
  for j in range(0,len(CDSstring)):
    CDSrow = CDSstring[j].strip('\n').split('\t')
    gene = CDSrow[5]
    start = int(CDSrow[1])
    stop = int(CDSrow[2])
    strand = CDSrow[3]
    if gene in genekeydict:
      CDSlist.append((gene,start,stop,strand))
    
  countvar=0
  
  #if wigstrand=='+':
  for line in CDSlist:
    if line[3]=='+':
      genename=line[0]
      start_index=line[1]
      stop_index=line[2]
      for index in range(start_index,stop_index+1):
        #if fnastring[index-2:index+2]=='TTGCG' and (fnastring[index-10:index-8]=='GG'):
        if (fnastring[index-3:index+2]=='TTGTG' or fnastring[index-3:index+2]=='TTGCG') and (fnastring[index-10:index-8]=='GG'):
          if index+1 in pause_dict:
            countvar+=1
            out_pausefile.write('>'+genename+'_'+line[3]+'_'+str(index+1)+'\n'+fnastring[index-39:index+2]+'\n')
            #out_pausefile.write(fnastring[index-35:index-5]+'\n')
          else:
            countvar+=1
            out_nopausefile.write('>'+genename+'_'+line[3]+'_'+str(index+1)+'\n'+fnastring[index-39:index+2]+'\n')
            #out_nopausefile.write(fnastring[index-35:index-5]+'\n')
            
  #if wigstrand=='-':
  #for line in CDSlist:
    elif line[3]=='-':
      genename=line[0]
      start_index=line[1]
      stop_index=line[2]
      for index in range(start_index,stop_index+1):
        #if (fnastring[index-1:index+4]=='CGCAA') and (fnastring[index+9:index+11]=='CC'):
        if (fnastring[index-1:index+4]=='CACAA' or fnastring[index-1:index+4]=='CGCAA') and (fnastring[index+9:index+11]=='CC'):
          pausestring=fnastring[index-1:index+40]
          pausestring=pausestring.replace('A','t').replace('C','g').replace('G','c').replace('T','a')
          pausestring=pausestring.replace('a','A').replace('c','C').replace('g','G').replace('t','T')
          if index+1 in pause_dict:
            countvar+=1
            out_pausefile.write('>'+genename+'_'+line[3]+'_'+str(index+1)+'\n'+pausestring[::-1]+'\n')
            #out_pausefile.write(pausestring[::-1]+'\n')
          else:
            countvar+=1
            out_nopausefile.write('>'+genename+'_'+line[3]+'_'+str(index+1)+'\n'+pausestring[::-1]+'\n')
            #out_nopausefile.write(pausestring[::-1]+'\n')
                        
  print countvar
  
  CDSfile.close()
  fnafile.close()
  pausefile.close()
  out_nopausefile.close()
  out_pausefile.close()
        
def main():
  ConsensusPauseFilter(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

if __name__ == '__main__':
  main()