import sys
import csv
import math
import numpy

def OccupancyAfterPause(wigfilename,genekeyfilename,pausefilename,outfilename,inputstrand):
  
  fnafile = open('c:\Python26\MG1655.fna', 'rU')
  fnastring=fnafile.read()
  fnastring=fnastring.replace('\n','')
  if inputstrand=='-':
    fnastring=fnastring.replace('A','t').replace('C','g').replace('G','c').replace('T','a')
    fnastring=fnastring.replace('a','A').replace('c','C').replace('g','G').replace('t','T')
  
  wigfile = open(wigfilename, 'rU') ##Use wigfile with all indices represented (e.g. greA_all.plus.wig)
  wigstring = wigfile.readlines()
  wigdict = {}
  
  genekeyfile=open(genekeyfilename)
  genekeydict={}
  
  pausefile=open(pausefilename, 'rU')
  pause_list=[]
  
  outfile = open(outfilename, 'wb')
  
  j=0
  for line in pausefile:
    if line[0]=='>':
      pause_row = line.strip('\n').strip('>').split('_')
      pause_index = int(pause_row[4])
      pause_fold = float(pause_row[-1])
      gene=pause_row[0]
      gene_avg = float(pause_row[-2])
      strand=pause_row[1]
      pause_list.append((pause_index,gene,strand,gene_avg))
      j+=1
  
  for line in genekeyfile:
    fields=line.replace('\n','').split('\t')
    gene=fields[0]
    transcriptionlevel=fields[2]
    genekeydict[gene]=transcriptionlevel
  
  TSSfile=open('MG1655_TSSs.mochiview.txt','rU')
  geneindexlist=[]
  for i,line in enumerate(TSSfile):
    #if i%1000==0:
      #print i
    
    fields=line.replace('\n','').split('\t')
    strand=fields[3]
    start=int(fields[1])
    stop=int(fields[2])
    seqID=(fields[12])
    genelength=stop-start
    
    if seqID in genekeydict and strand==inputstrand:
      geneindexlist.append((seqID,strand,start,stop,genekeydict[seqID]))
  
  for row in geneindexlist:
    gene=row[0]
    strand=row[1]
    start=row[2]
    stop=row[3]
    for i in range(start-50,stop+50):
      wigdict[i]=float(wigstring[i].strip('\n').split(' ')[1])
  
  distfrompause={}
  for i in range(-50,51):
    distfrompause[i]=0
  
  pauses_analyzed=0
  for row in pause_list:
    pause_index=row[0]
    gene=row[1]
    strand=row[2]
    gene_avg=row[3]
    if strand==inputstrand and gene in genekeydict:
    #if strand==inputstrand and gene=='ycaR':
      #print gene,pause_index,wigdict[pause_index-1]
      pauses_analyzed+=1
      for j in range(-50,51):
        count=wigdict[pause_index-1-j]
        distfrompause[j]+=count/gene_avg    
  
  j=50
  for element in distfrompause:
    distfrompause[j]/=pauses_analyzed
    outfile.write(str(j)+'\t'+str(distfrompause[j])+'\n')
    j-=1
  
  print 'pauses analyzed:',pauses_analyzed
  
  fnafile.close()
  TSSfile.close()
  wigfile.close()
  genekeyfile.close()
  outfile.close()  
  
def main():
  OccupancyAfterPause(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])

if __name__ == '__main__':
  main()