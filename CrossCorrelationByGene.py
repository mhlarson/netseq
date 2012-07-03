import sys
import csv
import numpy
import math

from numpy import *
# Runs a cross-correlation on two wig files in which each index is represented, even if there are zero counts at a given index.
# Riboseqfilename is moved forward/downstream across Netseqfilename over a given range. In bacteria, a transcriptional pause
# at a given sequence should lead to a ribosomal pause upstream of the transcriptional pause site. This program determines the 
# downstream shift in riboseqfilename that gives the maximal correlation. For plus strand, use as directed. For minus strand, reverse
# order of the input files. 

def CrossCorrelationByGene(ribowigfilename,netwigfilename,outcorrfilename,wigstrand):
  CDSfile = open ('CDS_MG1655.mochiview.txt','rU')
  CDSstring = CDSfile.readlines()
  CDSlist=[]
  
  ribowigfile = open(ribowigfilename, 'rU')
  netwigfile = open(netwigfilename,'rU')
  outcorrfile = open(outcorrfilename, 'wb')
  
  ribowigstring = ribowigfile.readlines()
  netwigstring = netwigfile.readlines()
  
  correlationgenome=[]
  
  for j in range(0,len(CDSstring)):
    CDSrow = CDSstring[j].strip('\n').split('\t')
    gene = CDSrow[5]
    start = int(CDSrow[1])
    stop = int(CDSrow[2])
    strand = CDSrow[3]
    CDSlist.append((gene,start,stop,strand))
        
  if wigstrand=='plus':
    for line in CDSlist:
      if line[3]=='+' and (line[2]-line[1]) > 50:
        genename=line[0]
        #print line[0]
        ribogene=[]
        netgene=[]
                
        for index in range(line[1]-1,line[2]):
          #print line[0],
          ribowigrow = ribowigstring[index].strip('\n').split(' ')
          ribowigcounts = float(ribowigrow[1])
          ribogene.append(ribowigcounts) 
          
          netwigrow = netwigstring[index].strip('\n').split(' ')
          netwigcounts = float(netwigrow[1])
          netgene.append(netwigcounts)
        
        if numpy.mean(ribogene)>0 and numpy.mean(netgene) >0:
          ribogene = (ribogene - numpy.mean(ribogene)) / (numpy.std(ribogene) * len(ribogene))
          netgene = (netgene - numpy.mean(netgene)) / numpy.std(netgene)
        
          ribogene=ribogene.tolist()
          netgene=netgene.tolist()
        
          correlationgene=[]
          for r in range(50):
            correlation = numpy.correlate(ribogene,netgene)
            ribogene.pop(-1)
            #x = ribogene.pop(-1)
            #ribogene.insert(0,x)
            correlationgene.append(correlation[0])
            r+=1
        
          correlationgenome.append(correlationgene)
          
  if wigstrand=='minus':
    for line in CDSlist:
      if line[3]=='-' and (line[2]-line[1]) > 50:
        genename=line[0]
        #print line[0]
        ribogene=[]
        netgene=[]
                
        for index in range(line[1]-1,line[2]):
          #print line[0],
          ribowigrow = ribowigstring[index].strip('\n').split(' ')
          ribowigcounts = float(ribowigrow[1])
          ribogene.append(ribowigcounts) 
          
          netwigrow = netwigstring[index].strip('\n').split(' ')
          netwigcounts = float(netwigrow[1])
          netgene.append(netwigcounts)
        
        if numpy.mean(ribogene)>0 and numpy.mean(netgene) >0:
          ribogene = (ribogene - numpy.mean(ribogene)) / (numpy.std(ribogene) * len(ribogene))
          netgene = (netgene - numpy.mean(netgene)) / numpy.std(netgene)
        
          ribogene=ribogene.tolist()
          netgene=netgene.tolist()
        
          correlationgene=[]
          for r in range(50):
            correlation = numpy.correlate(ribogene,netgene)
            ribogene.pop(-1)
            #x = ribogene.pop(-1)
            #ribogene.insert(0,0)
            correlationgene.append(correlation[0])
            r+=1
        
          correlationgenome.append(correlationgene)
  
  correlationarray = array(correlationgenome)
  correlationsum = correlationarray.sum(axis=0)/len(correlationarray)
  #correlationlist=correlationsum.tolist()
  
  correlationlist=[]
  for y in range(50):
    correlationlist.append((y,correlationsum[y]))
  
  print correlationlist
  
  #outcorrfile.write(correlationlist)
  
  csvlist = csv.writer(outcorrfile,delimiter = '\t')
  for row in correlationlist:
    csvlist.writerow(row)
  
  ribowigfile.close()
  netwigfile.close()
  outcorrfile.close()
        
def main():
  CrossCorrelationByGene(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

if __name__ == '__main__':
  main()