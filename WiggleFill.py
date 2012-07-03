import sys
import csv
import numpy
import math

def wigglefill(wigfilename,outwigfilename):
  wigfile = open(wigfilename, 'rU')
  wigstring = wigfile.readlines()
  wiglist = []
  outwigfile = open(outwigfilename, 'wb')
  
  for i in range(4639675):
    wiglist.append((i+1,0))
  
  for j in range(len(wigstring)-1):
    wigrow = wigstring[j].strip('\n').split(' ')
    index = int(wigrow[0])
    count = float(wigrow[1])
    wiglist[index-1]=(index,count)

  csvlist = csv.writer(outwigfile , delimiter = ' ')
  for row in wiglist:
    csvlist.writerow(row)
  
  wigfile.close()
  outwigfile.close()

def main():
  wigglefill(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()