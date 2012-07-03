import sys
import csv
import math
import numpy

def NextNucleotideAnalysis(wigfilename,genekeyfilename,inputstrand):
  
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
  
  A_energy=[]
  C_energy=[]
  G_energy=[]
  T_energy=[]
  
  #i=0
  for line in genekeyfile:
    #if i<336:
    fields=line.replace('\n','').split('\t')
    gene=fields[0]
    transcriptionlevel=fields[2]
    genekeydict[gene]=transcriptionlevel
    #i+=1
  
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
    for i in range(start,stop+1):
      wigdict[i]=float(wigstring[i].strip('\n').split(' ')[1])
  
  total_sequences=0
  
  #for j in range(1,21):
  for j in range(-20,21):  
    if j!=0:
      
      j_A=0
      j_C=0
      j_G=0
      j_T=0
      
      j_A_sem=0
      j_C_sem=0
      j_G_sem=0
      j_T_sem=0
      
      j_A_list=[]
      j_C_list=[]
      j_G_list=[]
      j_T_list=[]
  
      A_count=0
      C_count=0
      G_count=0
      T_count=0
  
      for row in geneindexlist:
        strand=row[1]
        gene=row[0]
        start=row[2]
        stop=row[3]
        for index in range(start+21,stop+1):
          incomingbase=fnastring[index-1]
          #if fnastring[index-12:index-10]=='GG' and (fnastring[index-5:index]=='TTGCG' or fnastring[index-5:index]=='TTGTG'):
          if incomingbase=='G':
            j_position=fnastring[index-1-j]
            count=wigdict[index-2]/float(genekeydict[gene])  ## Looking at index-1 for incoming base identity, but peak occurs at index-2
            total_sequences+=1
            if j_position=='A':
              A_count+=1
              j_A+=math.log(count+1)
              j_A_list.append(math.log(count+1))
            if j_position=='C':
              C_count+=1
              j_C+=math.log(count+1)
              j_C_list.append(math.log(count+1))
            if j_position=='G':
              G_count+=1
              j_G+=math.log(count+1)
              j_G_list.append(math.log(count+1))
            if j_position=='T':
              T_count+=1
              j_T+=math.log(count+1)
              j_T_list.append(math.log(count+1))
        
      j_A_sem=numpy.std(j_A_list)/math.sqrt(A_count)
      j_C_sem=numpy.std(j_C_list)/math.sqrt(C_count)
      j_G_sem=numpy.std(j_G_list)/math.sqrt(G_count)
      j_T_sem=numpy.std(j_T_list)/math.sqrt(T_count)
      
      j_A/=(A_count)
      j_C/=(C_count)
      j_G/=(G_count)
      j_T/=(T_count)
      
      j_N_mean=(j_A+j_C+j_G+j_T)/4
      
      j_A-=j_N_mean
      j_C-=j_N_mean
      j_G-=j_N_mean
      j_T-=j_N_mean
      
      A_energy.append((j_A,j_A_sem))
      C_energy.append((j_C,j_C_sem))
      G_energy.append((j_G,j_G_sem))
      T_energy.append((j_T,j_T_sem))
      
  outAfile = open('G_A.csv', 'wb')
  outCfile = open('G_C.csv', 'wb')
  outGfile = open('G_G.csv', 'wb')
  outTfile = open('G_T.csv', 'wb')
  
  for element in A_energy:
    outAfile.write(str(element[0])+'\n')
  for element in C_energy:
    outCfile.write(str(element[0])+'\n')
  for element in G_energy:
    outGfile.write(str(element[0])+'\n')
  for element in T_energy:
    outTfile.write(str(element[0])+'\n')
  
  total_sequences/=40
  print 'total sequences analyzed:',total_sequences
  
  fnafile.close()
  TSSfile.close()
  wigfile.close()
  genekeyfile.close()
  outAfile.close()
  outCfile.close()
  outGfile.close()
  outTfile.close()
              
def main():
  NextNucleotideAnalysis(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()