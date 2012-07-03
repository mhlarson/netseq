import sys
import csv
import numpy
import math

def pausefinder(netfilename):
  netfile = open(netfilename, 'rU')
    
  netstring = netfile.readlines()
  
  netlist=[]
    
  for i in range(len(netstring)):
    netrow = netstring[i].strip('\n').split(' ')
    netcounts = float(netrow[1])
    netlist.append(netcounts)  
  
  print numpy.mean(netlist)
  print numpy.std(netlist)
  
  #nethisto = numpy.histogram(netlist, bins=[0,10,20,30,40,50,60,70,80,90,100], density=False)
  #print nethisto[0][1]
  
#  csvlist = csv.writer(outwigfile , delimiter = ' ')
#  for row in wiglist:
#    csvlist.writerow(row)
   
  netfile.close()
  
def main():
  pausefinder(sys.argv[1])

if __name__ == '__main__':
  main()