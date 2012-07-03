import sys
import csv
import numpy
import math

# Runs a cross-correlation on two wig files in which each index is represented, even if there are zero counts at a given index.
# Riboseqfilename is moved forward/downstream across Netseqfilename over a given range. In bacteria, a transcriptional pause
# at a given sequence should lead to a ribosomal pause upstream of the transcriptional pause site. This program determines the 
# downstream shift in riboseqfilename that gives the maximal correlation. For plus strand, use as directed. For minus strand, reverse
# order of the input files. 

def correlate(riboseqfilename,netseqfilename,outcorrfilename):
  riboseqfile = open(riboseqfilename, 'rU')
  netseqfile = open(netseqfilename,'rU')
  
  outcorrfile = open(outcorrfilename, 'wb')
  
  riboseqstring = riboseqfile.readlines()
  netseqstring = netseqfile.readlines()

  riboseqlist=[]
  netseqlist=[]
  
  for i in range(len(riboseqstring)):
    riboseqrow = riboseqstring[i].strip('\n').split(' ')
    riboseqcounts = float(riboseqrow[1])
    riboseqlist.append(riboseqcounts)  
  for j in range(len(netseqstring)):  
    netseqrow = netseqstring[j].strip('\n').split(' ')
    netseqcounts = float(netseqrow[1])
    netseqlist.append(netseqcounts)
  
  #The next two lines normalize the input arrays prior to the cross-correlation.
  riboseqlist = (riboseqlist - numpy.mean(riboseqlist)) / (numpy.std(riboseqlist) * len(riboseqlist))
  netseqlist = (netseqlist - numpy.mean(netseqlist)) / numpy.std(netseqlist)
  
  riboseqlist=riboseqlist.tolist()
  netseqlist=netseqlist.tolist()
  
  correlationlist=[]
    
  for r in range(50):
    correlation = numpy.correlate(riboseqlist,netseqlist)
    x = riboseqlist.pop(-1)
    riboseqlist.insert(0,x)
    correlationlist.append((r,correlation[0]))
    r+=1
    
  maxcorr = max(correlationlist)
  
  csvlist = csv.writer(outcorrfile,delimiter = '\t')
  for row in correlationlist:
    csvlist.writerow(row)
  
  #print correlationlist
  #print 'Optimal Shift:',correlationlist.index(maxcorr)
  
  riboseqfile.close()
  netseqfile.close()

def main():
  correlate(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()