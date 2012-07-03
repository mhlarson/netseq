import sys
import os

def metagenelevels(inPlusFilename,inMinusFilename):
  
  outstartgenefile = open('metagenestartlevels.csv','wb')
  outendgenefile = open('metageneendlevels.csv','wb')
  
  GeneStartPlusReads={}
  GeneStartMinusReads={}
  GeneStartTotalReads={}
  
  GeneEndPlusReads={}
  GeneEndMinusReads={}
  GeneEndTotalReads={}
  
  for i in range(0,501):
    GeneStartPlusReads[i]=0
    GeneStartMinusReads[i]=0
    GeneStartTotalReads[i]=0
    
    GeneEndPlusReads[i]=0
    GeneEndMinusReads[i]=0
    GeneEndTotalReads[i]=0
  
  for i,line in enumerate(open('/home/mattlarson/genomes/MG1655_TSSs.mochiview.txt','rU')):
    if i%1000==0:
      print i,
        
    fields=line.replace('\n','').split('\t')
    strand=fields[3]
    start=int(fields[1])
    stop=int(fields[2])
    seqID=(fields[12])
    
    if stop-start>49: # and seqID!='yjtD':
      if strand=='+':
        inPlusFile=open(inPlusFilename, 'rU')
        filelist = inPlusFile.readlines()
        j=0
        for i in range(start,start+500,1):
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          GeneStartPlusReads[i]=reads
          j+=1
        j=0
        for i in range(stop,stop-500,-1)
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          GeneEndPlusReads[i]=reads
          j+=1
        
        #genetx=sum(geneReads[99:500])/400
        #geneavg=sum(geneReads)/len(geneReads)
        #outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
            
      elif strand=='-':
        inMinusFile=open(inMinusFilename, 'rU')
        filelist = inMinusFile.readlines()
        j=0
        for i in range(stop,stop-500,-1):
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          GeneStartMinusReads[j]+=reads
          j+=1
        j=0
        for in in range(start,start+500,1)
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          GeneEndMinusReads[i]=reads
          j+=1
      
        #genetx=sum(geneReads[99:500])/400
        #geneavg=sum(geneReads)/len(geneReads)
        #outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
                
  for i in range(0,500):
    GeneStartTotalReads[i]=GeneStartPlusReads[i]+GeneStartMinusReads[i] 
    GeneEndTotalReads[i]=GeneEndPlusReads[i]+GeneEndMinusReads[i] 
    outstartgenefile.write(str(i+1)+'\t'+str(GeneStartTotalReads[i])+'\n')
    outendgenefile.write(str(i+1)+'\t'+str(GeneEndTotalReads[i])+'\n')
  
  outstartgenefile.close()
  outendgenefile.close()


def main():
  metagenelevels(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()