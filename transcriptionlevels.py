import sys
import numpy as np
import os
import math

def transcriptionlevels(inPlusFilename,inMinusFilename):
  
  genome = {}
  for line in open('NC_000913_pause.fna', 'rU'):
    if line[0]=='>':
      chromosome=line[1:].rstrip('\r\n')
      genome[chromosome] = ''
    else:
      genome[chromosome]+=line.rstrip('\r\n')
      
  #numGenes=0
  outgenefile = open('txnlevels.csv','wb')
  
  for i,line in enumerate(open('MG1655_TSSs.mochiview.txt','rU')):
    if i%1000==0:
      print i
        
    fields=line.replace('\n','').split('\t')
    strand=fields[3]
    start=int(fields[1])
    stop=int(fields[2])
    seqID=(fields[12])
    if stop-start>49: # and seqID!='yjtD':
      if strand=='+':
        geneReads = []
        inFile=open(inPlusFilename, 'rU')
        filelist = inFile.readlines()
        stop_tolerance=0
        for i in range(start,stop+stop_tolerance):
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          geneReads.append(reads)
          
        genetx=sum(geneReads[99:500])/400
        geneavg=sum(geneReads)/len(geneReads)
        outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
            
      elif strand=='-':
        geneReads = []
        inFile=open(inMinusFilename, 'rU')
        filelist = inFile.readlines()
        stop_tolerance=0
        for i in range(start-stop_tolerance,stop):
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          geneReads.append(reads)
          
        genetx=sum(geneReads[99:500])/400
        geneavg=sum(geneReads)/len(geneReads)
        outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
                
                
  #print numGenes
  outgenefile.close()

def main():
  transcriptionlevels(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()