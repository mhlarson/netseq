import sys
import os

def metagenelevels(inPlusFilename,inMinusFilename,genekeyfilename):
  
  outstartgenefile = open('metagenestartlevels.csv','wb')
  outendgenefile = open('metageneendlevels.csv','wb')
   
  numgenes=0
  
  genekeyfile=open(genekeyfilename)
  genekeylist=[]
    
  for line in genekeyfile:
    fields=line.replace('\n','').split('\t')
    gene=fields[0]
    genekeylist.append(gene)
  
  GeneStartPlusReads={}
  GeneStartMinusReads={}
  GeneStartTotalReads={}
  
  GeneEndPlusReads={}
  GeneEndMinusReads={}
  GeneEndTotalReads={}
  
  for i in range(0,1001):
    GeneStartPlusReads[i]=0
    GeneStartMinusReads[i]=0
    GeneStartTotalReads[i]=0
    
    GeneEndPlusReads[i]=0
    GeneEndMinusReads[i]=0
    GeneEndTotalReads[i]=0
  
  for i,line in enumerate(open('MG1655_TSSs.mochiview.txt','rU')):
    if i%1000==0:
      print i
        
    fields=line.replace('\n','').split('\t')
    strand=fields[3]
    start=int(fields[1])
    stop=int(fields[2])
    seqID=(fields[12])
    genelength=stop-start
        
    if seqID in genekeylist and genelength>1000: # and seqID!='yjtD':
      #print seqID,
      if strand=='+':
        inPlusFile=open(inPlusFilename, 'rU')
        filelist = inPlusFile.readlines()
        geneReads = []
        
        j=0
        for i in range(start,start+500,1):
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          geneReads.append(reads)
          j+=1
          
        genetx=sum(geneReads[99:500])/400
        if genetx>1:
          numgenes+=1
          j=0
          for i in range(start,start+1000,1):
            stringline = filelist[i-1]
            countstring=stringline.split(' ')[1].strip('\n')
            if countstring =='': break
            reads = float(countstring) 
            GeneStartPlusReads[j]+=reads/genetx
            j+=1
                
          j=0
          for i in range(stop,stop-1000,-1):
            stringline = filelist[i-1]
            countstring=stringline.split(' ')[1].strip('\n')
            if countstring =='': break
            reads = float(countstring) 
            GeneEndPlusReads[j]+=reads/genetx
            j+=1
        
        #geneavg=sum(geneReads)/len(geneReads)
        #outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
            
      elif strand=='-':
        inMinusFile=open(inMinusFilename, 'rU')
        filelist = inMinusFile.readlines()
        geneReads = []
        
        j=0
        for i in range(stop,stop-500,-1):
          stringline = filelist[i-1]
          countstring=stringline.split(' ')[1].strip('\n')
          if countstring =='': break
          reads = float(countstring) 
          geneReads.append(reads)
          j+=1
        
        genetx=sum(geneReads[99:500])/400
        if genetx>1:
          j=0
          for i in range(stop,stop-1000,-1):
            stringline = filelist[i-1]
            countstring=stringline.split(' ')[1].strip('\n')
            if countstring =='': break
            reads = float(countstring) 
            geneReads.append(reads)
            GeneStartMinusReads[j]+=reads/genetx
            j+=1
        
          j=0
          for i in range(start,start+1000,1):
            stringline = filelist[i-1]
            countstring=stringline.split(' ')[1].strip('\n')
            if countstring =='': break
            reads = float(countstring) 
            GeneEndMinusReads[j]+=reads/genetx
            j+=1
        
        #geneavg=sum(geneReads)/len(geneReads)
        #outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
                
  for i in range(0,1000):
    GeneStartTotalReads[i]=(GeneStartPlusReads[i]+GeneStartMinusReads[i])/numgenes
    GeneEndTotalReads[i]=(GeneEndPlusReads[i]+GeneEndMinusReads[i])/numgenes
    outstartgenefile.write(str(i+1)+'\t'+str(GeneStartTotalReads[i])+'\n')
    outendgenefile.write(str(i+1)+'\t'+str(GeneEndTotalReads[i])+'\n')
  
  outstartgenefile.close()
  outendgenefile.close()

  print numgenes

def main():
  metagenelevels(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()