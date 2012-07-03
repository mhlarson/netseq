import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.stats as stats

def getPeakSeqs(genomeFilename,inPlusFilename,inMinusFilename,TSSIndex,seqLength,cutoff=0,peakStrengthMin=0,window=200,fold=4,peakmin=3,tnorm_thresh=2,rounds=10):
    pauseSeq=[]
    peakInfo=[]
    numGenes=0
    #load genome
    genome = loadGenome(genomeFilename)
    inPlus=open(inPlusFilename, 'r')
    inMinus=open(inMinusFilename,'r')
    for i,line in enumerate(open(TSSIndex)):
        if i%1000==0:
            print i
        
        fields=line.replace('\n','').split('\t')
        strand=fields[4]
        
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>500:
            if strand=='+':
                startIndex=int(fields[5])
                geneReads=getGeneReads(inPlus,startIndex,start,stop)
                
                genetx=sum(geneReads[99:500])/400
                if genetx >cutoff:
                    numGenes+=1
                    peaks,peakStrength = getGenePeaks(geneReads,window=window,fold=fold,peakmin=peakmin,tnorm_thresh=tnorm_thresh,rounds=rounds)
                    chrSeq=genome[fields[1]]
                    for i,peak in enumerate(peaks):
                        if peakStrength[i]>peakStrengthMin:
                            pauseSeq.append(findFlankSeq(chrSeq,start,peak[0],strand,seqLength))
                            peakInfo.append(fields+peak+[peakStrength[i],genetx])
            else:
                startIndex=int(fields[6])
                geneReads=getGeneReads(inMinus,startIndex,start,stop)
                
                genetx=sum(geneReads[-500:-99])/400
                if genetx >cutoff:
                    numGenes+=1
                    peaks,peakStrength = getGenePeaks(geneReads,window=window,fold=fold,peakmin=peakmin,tnorm_thresh=tnorm_thresh,rounds=rounds)
                    chrSeq=genome[fields[1]]
                    for i,peak in enumerate(peaks):
                        if peakStrength[i]>peakStrengthMin:
                            pauseSeq.append(findFlankSeq(chrSeq,start,peak[0],strand,seqLength))
                            peakInfo.append(fields+peak+[peakStrength[i],genetx])
    print numGenes
    return pauseSeq,peakInfo