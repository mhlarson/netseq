import numpy as np
#import matplotlib.pyplot as plt
import time
import scipy.stats as stats
#############################
# NET-seq analysis module
# Author: Stirling Churchman (modified by Matt Larson)
# Date: Decemeber 29, 2009
# Updated: May 4, 2011
############################

def getGeneReads(inFilename,start,stop,strand):
    offset=0
    geneReads = []
    inFile=open(inFilename, 'rU')
    filelist = inFile.readlines()
    stop_tolerance=0
    if strand=='+':
      for i in range(start,stop+stop_tolerance):
        stringline = filelist[i-1]
        countstring=stringline.split(' ')[1].strip('\n')
        if countstring =='': break
        reads = float(countstring) 
        geneReads.append(reads)
    elif strand=='-':
      for i in range(start-stop_tolerance,stop):
        stringline = filelist[i-1]
        countstring=stringline.split(' ')[1].strip('\n')
        if countstring =='': break
        reads = float(countstring) 
        geneReads.append(reads)
              
    return geneReads
    inFile.close()

def getChromReads(inFile,startIndex,start,stop):
    chromReads = []
    inFile.seek(startIndex)
    for i in range(start-1):line = inFile.readline()
    for i in range(start,stop+1):
        line = inFile.readline().replace('\n', '')
        reads = float(line)
        chromReads.append(reads)
    return chromReads

def getCorrelations(A,B):
    APlus=open(A+'_plus.txt', 'r')
    AMinus=open(A+'_minus.txt', 'r')
    BPlus=open(B+'_plus.txt', 'r')
    BMinus=open(B+'_minus.txt', 'r')
    BFields={}
    corr=[]
    Areads=[]
    Breads=[]
    geneNames=[]
    for line in open(B+'_index.txt'):
        fields = line.replace('\n','').split('\t')
        BFields[fields[0]]=fields
    for i,line in enumerate(open(A+'_index.txt')):
        if i%1000==0:
            print i
        #if i==500:break
        fields=line.replace('\n','').split('\t')
        strand=fields[4]
        
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>50:
            if strand=='+':
                startIndex=int(fields[5])
                BstartIndex=int(BFields[fields[0]][5])
                AGR=getGeneReads(APlus,startIndex,start,stop)
                BGR=getGeneReads(BPlus,BstartIndex,start,stop)
            else:
                startIndex=int(fields[6])
                BstartIndex=int(BFields[fields[0]][6])
                AGR=getGeneReads(AMinus,startIndex,start,stop)
                BGR=getGeneReads(BMinus,BstartIndex,start,stop)
            correlation = np.corrcoef(AGR,BGR)
            corr.append(correlation[0,1])
            Areads.append(sum(AGR[-500:])/500)
            Breads.append(sum(BGR[-500:])/500)
            geneNames.append(fields[0])
    return corr,Areads,Breads,geneNames

def getSpearCorrelations(A,B):
    APlus=open(A+'_plus.txt', 'r')
    AMinus=open(A+'_minus.txt', 'r')
    BPlus=open(B+'_plus.txt', 'r')
    BMinus=open(B+'_minus.txt', 'r')
    BFields={}
    corr=[]
    Areads=[]
    Breads=[]
    geneNames=[]
    for line in open(B+'_index.txt'):
        fields = line.replace('\n','').split('\t')
        BFields[fields[0]]=fields
    for i,line in enumerate(open(A+'_index.txt')):
        if i%1000==0:
            print i
        #if i==500:break
        fields=line.replace('\n','').split('\t')
        strand=fields[4]
        
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>500:
            if strand=='+':
                startIndex=int(fields[5])
                BstartIndex=int(BFields[fields[0]][5])
                AGR=getGeneReads(APlus,startIndex,start,stop)
                BGR=getGeneReads(BPlus,BstartIndex,start,stop)
            else:
                startIndex=int(fields[6])
                BstartIndex=int(BFields[fields[0]][6])
                AGR=getGeneReads(AMinus,startIndex,start,stop)
                BGR=getGeneReads(BMinus,BstartIndex,start,stop)
            correlation = stats.spearmanr(AGR,BGR)
            
            corr.append(correlation[0])
            Areads.append(sum(AGR[-500:])/500)
            Breads.append(sum(BGR[-500:])/500)
            geneNames.append(fields[0])
    return corr,Areads,Breads,geneNames

def loadGenome(genomeFilename):
    genome = {}
    for line in open(genomeFilename):
        if line[0]=='>':
            chromosome=line[1:].rstrip('\r\n')
            genome[chromosome] = ''
        else:
            genome[chromosome]+=line.rstrip('\r\n')
    return genome          
    
    
def getGenePeaks(GR,window,fold,peakmin,tnorm_thresh,rounds):
    #window=200,fold=4,peakmin=3,tnorm_thresh=2,rounds=10
    gene_reads = [num for num in GR]
    peaks = []
    peakStrength = []
    i=0
    numPeaks=-1
    while i<rounds and len(peaks)>numPeaks:
        numPeaks=len(peaks)
        #sliding window
        # remove peaks found in previous round and set them to the average value
        #for that region
        for peak in peaks:
            gene_reads[peak[0]]=tnorm[peak[0]]
            
        
        delta1 = []
        delta2 = []
        genelength = len(gene_reads)
        win_starts = [0]*(window/2)
        win_starts+=range( genelength-(window) )
        win_starts+=[genelength-window]*(window/2)
        win_stops=[window]*(window/2)
        win_stops += range((window),genelength)
        win_stops+=[genelength]*(window/2)
        
        tnorm = [np.mean([num for num in gene_reads[win_starts[index]:win_stops[index]] if num >0]) for index,peak in enumerate(gene_reads)]
        #tnorm = [np.median(gene_reads[win_starts[index]:win_stops[index]]) for index,peak in enumerate(gene_reads)]
      
        #thresh = [tnorm[index]+fold*np.std(gene_reads[win_starts[index]:win_stops[index]]) for index,peak in enumerate(gene_reads)]
        thresh =[]
        for index,peak in enumerate(gene_reads):
            geneReadsGTzero = [num for num in gene_reads[win_starts[index]:win_stops[index]] if num>0]
            thresh.append(tnorm[index]+fold*np.std(geneReadsGTzero))
        
        
        peaks+= [[index,peak] for (index,peak) in enumerate(gene_reads) if peak >thresh[index] and tnorm[index]>tnorm_thresh and peak>peakmin]
        #numPeaks = len(peaks)
        i+=1
    for peak in peaks:
        maxtnorm=max(tnorm)
        #peakStrength.append(peak[1]/maxtnorm)
        peakStrength.append(peak[1]/tnorm[peak[0]])
    #print i
    #genetx1=sum(GR[99:500])/400
    #genetx2=sum(GR[-500:-99])/400
    #print genetx1
    #print genetx2
    ##if len(peaks)>2:
    #if genetx1>=0 or genetx2 >=0:
    #    
    ##if i >0:
    #    plt.figure(1)
    #    plt.plot(GR)
    ###    #
    #    plt.plot([peak[0] for peak in peaks],[peak[1] for peak in peaks], 'g*')
    #    #plt.figure(2)
    #    #loggedGR=np.log(np.array(GR))
    #    #plt.hist(loggedGR)
    #    #GRA=np.array(GR)
    #    #sp = np.fft.fft(np.array(GR))
    #    #freq = np.fft.fftfreq(GRA.shape[-1])
    #    #plt.plot(freq,sp)
    #    #plt.xlim(0,.4)
    #    plt.show()    
    return peaks,peakStrength

def getPeakSeqs(genomeFilename,inPlusFilename,inMinusFilename,TSSIndex,seqLength,cutoff,peakStrengthMin,window,fold,peakmin,tnorm_thresh,rounds,promoterfile):
    
    outgenefile = open('txnlevels.csv','wb')
    
    pauseSeq=[]
    peakInfo=[]
    promdict={}
    numGenes=0
    #load genome
    genome = loadGenome(genomeFilename)
    #inPlus=open(inPlusFilename, 'rU')
    #inMinus=open(inMinusFilename,'rU')
    
    for line in (open(promoterfile)):
      promfields=line.rstrip('\r\n').split('\t')
      promseqID=promfields[12].strip('1').strip('2').strip('3').strip('4').strip('5').strip('6').strip('7')[:-1]
      promstart=promfields[1]
      promstrand=promfields[3]
      
      if promseqID in promdict:
        if promstrand=="+":
          if promstart<promdict[promseqID]:
            promdict[promseqID]=promstart
        
        if promstrand=="-":
          if promstart>promdict[promseqID]:
            promdict[promseqID]=promstart
          
      else:
        promdict[promseqID]=promstart
      
    for i,line in enumerate(open(TSSIndex)):
        if i%1000==0:
            print i
        
        fields=line.replace('\n','').split('\t')
        strand=fields[3]
        start=int(fields[1])
        stop=int(fields[2])
        seqID=(fields[12])
        if stop-start>49: # and seqID!='yjtD':
            if strand=='+':
                #if seqID in promdict:
                #  start=int(promdict[seqID])
                  
                geneReads=getGeneReads(inPlusFilename,start,stop,strand)
                                
                genetx=sum(geneReads[99:500])/400
                geneavg=sum(geneReads)/len(geneReads)
                outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
                if genetx >cutoff:
                    numGenes+=1
                    peaks,peakStrength = getGenePeaks(geneReads,window=window,fold=fold,peakmin=peakmin,tnorm_thresh=tnorm_thresh,rounds=rounds)
                    chrSeq=genome['MG1655']
                    startlist=[]
                    startlist.append(start)
                    for i,peak in enumerate(peaks):
                        if peakStrength[i]>peakStrengthMin:
                            pauseSeq.append(findFlankSeq(chrSeq,start,peak[0],strand,seqLength))
                            peakInfo.append(fields+peak+startlist+[peakStrength[i],geneavg])
            elif strand=='-':
                #if seqID in promdict:
                #  stop=int(promdict[seqID])
                  
                geneReads=getGeneReads(inMinusFilename,start,stop,strand)
                
                genetx=sum(geneReads[-500:-99])/400
                geneavg=sum(geneReads)/len(geneReads)
                outgenefile.write(seqID+'_'+strand+'_'+str(geneavg)+'\n')
                if genetx >cutoff:
                    numGenes+=1
                    peaks,peakStrength = getGenePeaks(geneReads,window=window,fold=fold,peakmin=peakmin,tnorm_thresh=tnorm_thresh,rounds=rounds)
                    chrSeq=genome['MG1655']
                    stoplist=[]
                    stoplist.append(stop)
                    for i,peak in enumerate(peaks):
                        if peakStrength[i]>peakStrengthMin:
                            pauseSeq.append(findFlankSeq(chrSeq,start,peak[0],strand,seqLength))
                            peakInfo.append(fields+peak+stoplist+[peakStrength[i],geneavg])
    print numGenes
    #print 'pauseSeq',pauseSeq
    #print 'peakInfo',peakInfo
    return pauseSeq,peakInfo

def getPeakTime(inPlusFilename,inMinusFilename,TSSIndex,seqLength,cutoff=0,peakStrengthMin=0,window=200,fold=4,peakmin=3,tnorm_thresh=2,rounds=10):
    pauseFraction=[]
    numGenes=0
    
    inPlus=open(inPlusFilename, 'r')
    inMinus=open(inMinusFilename,'r')
    for i,line in enumerate(open(TSSIndex)):
        if i%1000==0:
            print i
        #if i>100:
        #    break
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
                    peakSum = sum([num[1] for num in peaks])
            
                    pauseFraction.append([fields[0],(peakSum/sum(geneReads)),genetx])
                    
            
            else:
                startIndex=int(fields[6])
                geneReads=getGeneReads(inMinus,startIndex,start,stop)
                genetx=sum(geneReads[-500:-99])/400
                if genetx >cutoff:
                    numGenes+=1
                    peaks,peakStrength = getGenePeaks(geneReads,window=window,fold=fold,peakmin=peakmin,tnorm_thresh=tnorm_thresh,rounds=rounds)
                    peaks,peakStrength = getGenePeaks(geneReads,window=window,fold=fold,peakmin=peakmin,tnorm_thresh=tnorm_thresh,rounds=rounds)
                    peakSum = sum([num[1] for num in peaks])
            
                    pauseFraction.append([fields[0],(peakSum/sum(geneReads)),genetx])
                    
            
    print numGenes
    return pauseFraction

def getUniquePeakSeqs(genomeFilename,WT,WTRep,DST1,seqLength=15,cutoff=0,peakStrengthMin=0,WTpeakmin=3,WTRepPeakMin=5,DST1PeakMin=5,rounds=10):
    pauseSeq=[]
    peakInfo=[]
    WTRepFields={}
    DST1Fields={}
    numGenes=0
    plot=0
    #load genome
    print 'Loading genome'
    genome = loadGenome(genomeFilename)
    print 'Opening data files'
    WTPlus=open(WT+'_plus.txt', 'r')
    WTMinus=open(WT+'_minus.txt', 'r')
    WTRepPlus=open(WTRep+'_plus.txt', 'r')
    WTRepMinus=open(WTRep+'_minus.txt', 'r')
    DST1Plus=open(DST1+'_plus.txt', 'r')
    DST1Minus=open(DST1+'_minus.txt', 'r')
    for line in open(WTRep+'_index.txt'):
        fields = line.replace('\n','').split('\t')
        WTRepFields[fields[0]]=fields
    for line in open(DST1+'_index.txt'):
        fields = line.replace('\n','').split('\t')
        DST1Fields[fields[0]]=fields
    for i,line in enumerate(open(WT+'_index.txt')):
        if i%1000==0:
            print i
        fields=line.replace('\n','').split('\t')
        strand=fields[4]
        
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>500:
            if strand=='+':
                startIndex=int(fields[5])
                WTRepstartIndex=int(WTRepFields[fields[0]][5])
                DST1startIndex=int(DST1Fields[fields[0]][5])
                WTGR=getGeneReads(WTPlus,startIndex,start,stop)
                WTRepGR=getGeneReads(WTRepPlus,WTRepstartIndex,start,stop)
                DST1GR=getGeneReads(DST1Plus,DST1startIndex,start,stop)
                
            else:
                startIndex=int(fields[6])
                WTRepstartIndex=int(WTRepFields[fields[0]][6])
                DST1startIndex=int(DST1Fields[fields[0]][6])
                WTGR=getGeneReads(WTMinus,startIndex,start,stop)
                WTRepGR=getGeneReads(WTRepMinus,WTRepstartIndex,start,stop)
                DST1GR=getGeneReads(DST1Minus,DST1startIndex,start,stop)
            reads = sum(WTGR[:100])/100
            if reads >int(cutoff):
                
                numGenes+=1
                WTpeaks,WTpeakStrength = getGenePeaks(WTGR,peakmin=WTpeakmin,rounds=rounds)
                WTReppeaks,WTReppeakStrength = getGenePeaks(WTRepGR,peakmin=WTRepPeakMin,rounds=rounds)
                DST1peaks,DST1peakStrength = getGenePeaks(DST1GR,peakmin=DST1PeakMin,rounds=rounds)
                WTpeakPos = [peak[0] for i,peak in enumerate(WTpeaks) if WTpeakStrength[i]>=int(peakStrengthMin)]
                WTReppeakPos = [peak[0] for peak in WTReppeaks]
                DST1peakPos = [peak[0] for peak in DST1peaks]
                #print len(DST1peakPos)
                #print len(WTpeakPos)
                #RepPeaks = [peak for i,peak in enumerate(WTpeaks) if (peak[0] in WTReppeakPos and peak[0] not in DST1peakPos)]
                RepPeaks = [peak for i,peak in enumerate(WTpeakPos) if (peak not in DST1peakPos)]

                chrSeq=genome[fields[1]]
                for i,peak in enumerate(RepPeaks):
                    #if peakStrength[i]>peakStrengthMin:
                    pauseSeq.append(findFlankSeq(chrSeq,start,peak,strand,seqLength))
                    peakInfo.append(fields+[peak])
                if plot:
                    fig = plt.figure(1)
                    ax = fig.add_subplot(111)
                    x=range(start,stop+1)
                    geneReads22=[-num for num in DST1GR]
                    plt.bar(x,WTGR,color='red',edgecolor='red')
                    plt.plot([peak[0]+start for peak in WTpeaks],[peak[1] for peak in WTpeaks], 'g*')
                    plt.bar(x,geneReads22,color='blue',edgecolor='blue') #DST1
                    plt.plot([peak[0]+start for peak in DST1peaks],[-peak[1] for peak in DST1peaks], 'g*')
                    plt.title(fields[0])
                    plt.grid()
                    if RepPeaks !=[]:
                        #print RepPeaks
                        plt.plot([i+start for i in RepPeaks],[0]*len(RepPeaks),'*k')
                    fig = plt.figure(2)
                    ax = fig.add_subplot(111)
                    x=range(start,stop+1)
                    geneReads22=[-num for num in WTRepGR]
                    plt.bar(x,WTGR,color='red',edgecolor='red')
                    plt.plot([peak[0]+start for peak in WTpeaks],[peak[1] for peak in WTpeaks], 'g*')
                    plt.bar(x,geneReads22,color='blue',edgecolor='blue') #DST1
                    plt.plot([peak[0]+start for peak in WTReppeaks],[-peak[1] for peak in WTReppeaks], 'g*')
                    plt.title(fields[0])
                    plt.grid()
                    plt.show()
    print numGenes
    return pauseSeq,peakInfo

                           
def reverseCompliment(seq):
    compliment = {'A':'T','T':'A','C':'G','G':'C'}
    seq = list(seq)
    seq.reverse()
    seqrc = ''
    for s in seq:
        seqrc = seqrc+compliment[s]
    return seqrc

def findFlankSeq(chr_seq,start,index,strand,window):
    if strand == '+':
        upSeq = chr_seq[start+index-window:start+index]
        downSeq = chr_seq[start+index:start+index+window]
    else:
        downSeq = reverseCompliment(chr_seq[start+index-window-1:start+index-1])
        upSeq = reverseCompliment(chr_seq[start+index-1:start+index+window-1])
    return upSeq,downSeq



def getSATandemSlide(inputFile1,tssFile,excludeGenes,winSense,winAnti,delta,win_slide_anti,win_slide_sense):
     # open the track files for the plus reads
    filename1 = inputFile1+'_plus.txt'
    filename2 = inputFile1+'_minus.txt'
    startIndex={}
    for line in open(inputFile1+'_index.txt'):
        fields=line.replace('\n','').split('\t')
        startIndex[fields[0]]=[int(fields[5]),int(fields[6])]

    inFilePlus = file(filename1,'r')
    inFileMinus = file(filename2,'r')
    senseReads=[]
    antiReads=[]
    geneNames=[]
    Adelta = []
    Sdelta = []
    ratio = []
    
    flags=[]
    all_genes = {}
    count = 0
    step = 10
    plot = 0
    for i,line in enumerate(open(tssFile)):
        pair_reads_plus = []
        pair_reads_minus = []
        pairinfo = line.replace('\n','').split('\t')
        if i%1000 ==0:
            print i
        if pairinfo[0]=='flag':
            pairinfo.pop(0)
            flag=1
        else:
            flag=0
        start = int(pairinfo[1])
        stop = int(pairinfo[7])
        gene1_name = pairinfo[4]
        gene1_length = int(pairinfo[2])-start
        gene2_length = stop-int(pairinfo[6])
        if stop-start>winSense+winAnti+delta+win_slide_sense+win_slide_anti and gene1_length>winSense+win_slide_sense and gene2_length>winSense+win_slide_sense:
            #print pairinfo
            
            pair_reads_plus = getGeneReads(inFilePlus,startIndex[gene1_name][0],start,stop)
            pair_reads_minus = getGeneReads(inFileMinus,startIndex[gene1_name][1],start,stop)
            
            current_chromosome = pairinfo[0]
            flags.append(flag)
        else: continue
        if pairinfo[3]=='+':
            if pairinfo[9] in excludeGenes:continue
            startp = int(pairinfo[6])-start
            
            sense = [sum(pair_reads_plus[sp:sp+winSense]) for sp in range(startp,startp+win_slide_sense,step)]
            anti = [sum(pair_reads_minus[startp-d-winAnti:startp-d]) for d in range(delta,delta+win_slide_anti,step)]
            I = np.argsort(np.array(sense))
            IA = np.argsort(np.array(anti))
            Adelta.append(delta+IA[-1]*step)
            Sdelta.append(I[-1]*step)
            
            sense = max(sense)
            anti = max(anti)
            
            if float(sense)>0:
                ratio.append(float(anti)/sense)
            else:
                ratio.append(-1)
            senseReads.append(sense)
            antiReads.append(anti)
            geneNames.append(pairinfo[9])
            if plot:
                print I[-1]
                print sense
                print anti
                print pairinfo[9]
                plt.figure()
                plt.plot(pair_reads_plus,'b',[-reads for reads in pair_reads_minus], 'r')
                plt.vlines(int(pairinfo[6])-start,0,max(pair_reads_plus),linestyle='dashed')
                plt.vlines(int(pairinfo[6])-start+(I[-1]*step),0,max(pair_reads_plus),linestyle='dashed')
                plt.vlines(int(pairinfo[6])-start+(I[-1]*step)+winSense,0,max(pair_reads_plus),linestyle='dashed')
                plt.vlines(int(pairinfo[6])-start-delta-(IA[-1]*step),0,-max(pair_reads_minus),linestyle='dashed')
                plt.vlines(int(pairinfo[6])-start-delta-(IA[-1]*step)-winAnti,0,-max(pair_reads_minus),linestyle='dashed')
                
                plt.title(pairinfo[9])
                plt.show()
        else:
            if pairinfo[4] in excludeGenes:continue
            stopp = int(pairinfo[2])-start
            anti=[sum(pair_reads_plus[stopp+d:stopp+d+winAnti]) for d in range(delta,delta+win_slide_anti,step)]
            sense = [sum(pair_reads_minus[sp-winSense:sp]) for sp in range(stopp-win_slide_sense,stopp,step)]
            I = np.argsort(np.array(sense))
            IA = np.argsort(np.array(anti))
            Adelta.append(delta+IA[-1]*step)
            Sdelta.append(I[-1]*step)
            sense = max(sense)
            anti = max(anti)
           
            if float(sense)>0:
                ratio.append(float(anti)/sense)
            else:
                ratio.append(-1)
            senseReads.append(sense)
            antiReads.append(anti)
            geneNames.append(pairinfo[4])
            if plot:
                print I[-1]
                print sense
                print anti
                print pairinfo[4]
                plt.figure()
                #pair_reads_plus.reverse()
                #pair_reads_minus.reverse()
                
                plt.plot(pair_reads_plus,'b',[-reads for reads in pair_reads_minus], 'r')
                plt.vlines(int(pairinfo[2])-start,0,-max(pair_reads_minus),linestyle='dashed')
                plt.vlines(int(pairinfo[2])-start-(step*I[-1]),0,-max(pair_reads_minus),linestyle='dashed')
                plt.vlines(int(pairinfo[2])-start-(step*I[-1])-winSense,0,-max(pair_reads_minus),linestyle='dashed')
                plt.vlines(int(pairinfo[2])-start+delta+(step*IA[-1]),0,max(pair_reads_plus),linestyle='dashed')
                plt.vlines(int(pairinfo[2])-start+delta+(step*IA[-1])+winAnti,0,max(pair_reads_plus),linestyle='dashed')
                plt.title(pairinfo[4])
                plt.show()
                
    return senseReads,antiReads,ratio,geneNames,flags,Sdelta,Adelta
def getSATandemSlideUp(inputFile1,tssFile,excludeGenes,winSense,winAnti,delta,win_slide):
     # open the track files for the plus reads
    filename1 = inputFile1+'_plus.txt'
    filename2 = inputFile1+'_minus.txt'
    startIndex={}
    for line in open(inputFile1+'_index.txt'):
        fields=line.replace('\n','').split('\t')
        startIndex[fields[0]]=[int(fields[5]),int(fields[6])]

        
    
    inFilePlus = file(filename1,'r')
    inFileMinus = file(filename2,'r')
    # readout the first line since it won't be needed
    
    senseReads=[]
    antiReads=[]
    geneNames=[]
    ratio = []
    win=[]
    win2=[]
    all_genes = {}
    count = 0
    for i,line in enumerate(open(tssFile)):
        pair_reads_plus = []
        pair_reads_minus = []
        pairinfo = line.replace('\n','').split('\t')
        if i%1000 ==0:
            print i
        
        start = int(pairinfo[1])
        stop = int(pairinfo[7])
        gene1_name = pairinfo[4]
        gene1_length = int(pairinfo[2])-start
        gene2_length = stop-int(pairinfo[6])
        if stop-start>winSense+winAnti+delta+win_slide and gene1_length>winSense+win_slide and gene2_length>winSense+win_slide:
            #print pairinfo
            
            pair_reads_plus = getGeneReads(inFilePlus,startIndex[gene1_name][0],start,stop)
            pair_reads_minus = getGeneReads(inFileMinus,startIndex[gene1_name][1],start,stop)
            
            current_chromosome = pairinfo[0]
        else: continue
        if pairinfo[3]=='+':
            if pairinfo[9] in excludeGenes:continue
            startp = int(pairinfo[6])-start
            startu=0
            sense = [sum(pair_reads_plus[sp:sp+winSense]) for sp in range(startu,startu+win_slide,10)]
            anti = [sum(pair_reads_minus[startp-d-winAnti:startp-d]) for d in range(delta,delta+win_slide,10)]
            I = np.argsort(np.array(sense))
            win.append(I[-1])
            
            sense = max(sense)
            #print sense
            anti = max(anti)
            
            if float(sense)>0:
                ratio.append(float(anti)/sense)
            else:
                ratio.append(-1)
            senseReads.append(sense)
            antiReads.append(anti)
            geneNames.append(pairinfo[9])
            #print I[-1]
            #print sense
            #print anti
            #print pairinfo[9]
            #plt.figure()
            #plt.plot(pair_reads_plus,'b',[-reads for reads in pair_reads_minus], 'r')
            #plt.title(pairinfo[9])
            #plt.show()
        else:
            if pairinfo[4] in excludeGenes:continue
            
            stopp = int(pairinfo[2])-start
            stopu=stop-start
            anti=[sum(pair_reads_plus[stopp+d:stopp+d+winAnti]) for d in range(delta,delta+win_slide,10)]
            sense = [sum(pair_reads_minus[sp-winSense:sp]) for sp in range(stopu-win_slide,stopu,10)]
            I = np.argsort(np.array(sense))
            win2.append(I[-1])
            
            sense = max(sense)
            #print sense
            anti = max(anti)
           
            if float(sense)>0:
                ratio.append(float(anti)/sense)
            else:
                ratio.append(-1)
            senseReads.append(sense)
            antiReads.append(anti)
            geneNames.append(pairinfo[4])
            #print I[-1]
            #print sense
            #print anti
            #print pairinfo[4]
            #plt.figure()
            #plt.plot(pair_reads_plus,'b',[-reads for reads in pair_reads_minus], 'r')
            #plt.title(pairinfo[4])
            #plt.show()
            #if anti ==0:
            #    print sense
            #    print anti
            #    print pairinfo
            #    plt.figure()
            #    plt.plot(pair_reads_plus,'b',[-reads for reads in pair_reads_minus], 'r')
            #    plt.title(pairinfo[4])
            #    plt.show()
            #if pairinfo[4]=='YDL193W':
            #    
            #    plt.figure()
            #    plt.plot(pair_reads_plus,'b',pair_reads_minus, 'r')
            #    plt.title(pairinfo[4])
            #    plt.show()
    return senseReads,antiReads,ratio,geneNames,win,win2
def getAverageGene(inFileBase,length,minGenelength,lower=0,upper = 100000, normalize=1,genes=[]):
    inFilePlus = open(inFileBase+'_plus.txt','r')
    inFileMinus = open(inFileBase+'_minus.txt','r')
    averageGenePos={}
    averageGene=[]
    geneCounts={}
    geneRejects=[]
    geneNum=0
    for i,line in enumerate(open(inFileBase+'_index.txt')):
        if i%1000==0:print i
        fields=line.replace('\n','').split('\t')
        if fields[0] in genes:
            geneRejects.append(fields[0])
            continue
        
        strand=fields[4]
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>=minGenelength:
            if strand == '+':
                startIndex=int(fields[5])
                geneReads = getGeneReads(inFilePlus,startIndex,start,stop)
            else:
                startIndex=int(fields[6])
                geneReads = getGeneReads(inFileMinus,startIndex,start,stop)
                geneReads.reverse()
            norm = float(sum(geneReads[-500:]))/500
        
            if norm>lower and norm<upper:
                if not normalize:norm=1
                geneNum+=1
                for pos,read in enumerate(geneReads):
                    
                    if pos in averageGenePos:
                        averageGenePos[pos]+=read/norm
                        geneCounts[pos]+=1
                    else:
                        averageGenePos[pos]=read/norm
                        geneCounts[pos] = 1
            else: geneRejects.append(fields[0])
        else: geneRejects.append(fields[0])
                        
            
    positions = averageGenePos.keys()
    positions.sort()
    for pos in positions:
        averageGene.append(float(averageGenePos[pos])/geneCounts[pos])
    print 'number of genes: %s' % (geneNum)
    return averageGene,geneRejects
def getAverageGeneComp(inFileBase1,inFileBase2,length,minGenelength,lower=0,upper = 100000, normalize=1,genes=[]):
    inFilePlus = open(inFileBase1+'_plus.txt','r')
    inFileMinus = open(inFileBase1+'_minus.txt','r')
    inFilePlus2 = open(inFileBase2+'_plus.txt','r')
    inFileMinus2 = open(inFileBase2+'_minus.txt','r')
    averageGenePos={}
    averageGene=[]
    averageGenePos2={}
    averageGene2=[]
    geneCounts={}
    geneIndex2={}
    geneNum=0
    for i,line in enumerate(open(inFileBase2+'_index.txt')):
        fields=line.replace('\n','').split('\t')
        geneIndex2[fields[0]]=[int(fields[5]),int(fields[6])]
    for i,line in enumerate(open(inFileBase1+'_index.txt')):
        if i%1000==0:print i
        fields=line.replace('\n','').split('\t')
        if fields[0] not in geneIndex2:continue
        if fields[0] in genes:continue
        
        strand=fields[4]
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>=minGenelength:
            
            if strand == '+':
                startIndex=int(fields[5])
                geneReads = getGeneReads(inFilePlus,startIndex,start,stop)
                startIndex2=geneIndex2[fields[0]][0]
                geneReads2 = getGeneReads(inFilePlus2,startIndex2,start,stop)
            else:
                startIndex=int(fields[6])
                geneReads = getGeneReads(inFileMinus,startIndex,start,stop)
                startIndex2=geneIndex2[fields[0]][1]
                geneReads2 = getGeneReads(inFileMinus2,startIndex2,start,stop)
                geneReads.reverse()
                geneReads2.reverse()
            
            norm = float(sum(geneReads[100:501]))/len(geneReads[100:500])
            norm2 = float(sum(geneReads2[100:501]))/len(geneReads[100:500])
            
            if (norm>lower and norm<upper) and (norm2>lower and norm2<upper):
                
                if not normalize:
                    norm=1
                    norm2=1
                geneNum+=1
                for pos,read in enumerate(geneReads):
                    read2 = geneReads2[pos]
                    if pos in averageGenePos:
                        averageGenePos[pos]+=read/norm
                        geneCounts[pos]+=1
                    else:
                        averageGenePos[pos]=read/norm
                        geneCounts[pos] = 1
                    if pos in averageGenePos2:
                        averageGenePos2[pos]+=read2/norm2
                        
                    else:
                        averageGenePos2[pos]=read2/norm2
                        
                   
    positions = averageGenePos.keys()
    positions.sort()
    for pos in positions:
        averageGene.append(float(averageGenePos[pos])/geneCounts[pos])
    positions = averageGenePos2.keys()
    positions.sort()
    for pos in positions:
        averageGene2.append(float(averageGenePos2[pos])/geneCounts[pos])
    print geneNum
    return averageGene,averageGene2
def getAverageTandemGene(inputFile1,tssFile,minGenelength,lower=0,upper = 100000, normalize=1):
     # open the track files for the plus reads
    filename1 = inputFile1+'_plus.txt'
    filename2 = inputFile1+'_minus.txt'
    startIndex={}
    for line in open(inputFile1+'_index.txt'):
        fields=line.replace('\n','').split('\t')
        startIndex[fields[0]]=[int(fields[5]),int(fields[6])]

        
    
    inFilePlus = file(filename1,'r')
    inFileMinus = file(filename2,'r')
    # readout the first line since it won't be needed
    
    senseReads=[]
    antiReads=[]
    geneNames=[]
    ratio = []
    win=[]
    win2=[]
    all_genes = {}
    count = 0
    averageGenePosPlus={}
    averageGenePosMinus={}
    averageGenePlus=[]
    averageGeneMinus=[]
    geneCounts={}
    geneRejects=[]
    geneNum=0
    for i,line in enumerate(open(tssFile)):
        pair_reads_plus = []
        pair_reads_minus = []
        pairinfo = line.replace('\n','').split('\t')
        if i%1000 ==0:
            print i
        if pairinfo[0]=='flag':
            pairinfo.pop(0)
        print pairinfo
        start = int(pairinfo[1])
        stop = int(pairinfo[7])
        gene1_name = pairinfo[4]
        gene1_length = int(pairinfo[2])-start
        gene2_length = stop-int(pairinfo[6])
        if gene1_length>minGenelength and gene2_length>minGenelength:
            #print pairinfo
            
            pair_reads_plus = getGeneReads(inFilePlus,startIndex[gene1_name][0],start,stop)
            pair_reads_minus = getGeneReads(inFileMinus,startIndex[gene1_name][1],start,stop)
            
            
        else: continue
        if pairinfo[3]=='-':
            startp=stop-int(pairinfo[2])
            pair_reads_plus.reverse()
            pair_reads_minus.reverse()
            temp=pair_reads_plus
            pair_reads_plus=pair_reads_minus
            pair_reads_minus=temp
        else:
            startp=int(pairinfo[6])-start
            
            
        norm = float(sum(pair_reads_minus[startp-800:startp-200]))/600
    
    
        if norm>lower and norm<upper:
                if not normalize:norm=1
                geneNum+=1
                for index,read in enumerate(pair_reads_minus):
                    readPlus=pair_reads_plus[index]
                    pos=index-startp
                    if pos in averageGenePosPlus:
                        averageGenePosMinus[pos]+=read/norm
                        averageGenePosPlus[pos]+=readPlus/norm
                        geneCounts[pos]+=1
                    else:
                        averageGenePosMinus[pos]=read/norm
                        averageGenePosPlus[pos]=readPlus/norm
                        geneCounts[pos] = 1
                
    positions = averageGenePosPlus.keys()
    positions.sort()
    for pos in positions:
        averageGenePlus.append(float(averageGenePosPlus[pos])/geneCounts[pos])
        averageGeneMinus.append(float(averageGenePosMinus[pos])/geneCounts[pos])
    #positions = averageGenePosMinus.keys()
    #positions.sort()
    
        
    print geneNum
    return positions,averageGenePlus,averageGeneMinus
def getAverageTandemGeneComp(inputFile1,inputFile2,tssFile,minGenelength,lower=0,upper = 100000, normalize=1):
     # open the track files for the plus reads
    filename1 = inputFile1+'_plus.txt'
    filename2 = inputFile1+'_minus.txt'
    filename3 = inputFile2+'_plus.txt'
    filename4 = inputFile2+'_minus.txt'
    startIndex={}
    for line in open(inputFile1+'_index.txt'):
        fields=line.replace('\n','').split('\t')
        startIndex[fields[0]]=[int(fields[5]),int(fields[6])]
    startIndex2={}
    for line in open(inputFile2+'_index.txt'):
        fields=line.replace('\n','').split('\t')
        startIndex2[fields[0]]=[int(fields[5]),int(fields[6])]
    
    inFilePlus = file(filename1,'r')
    inFileMinus = file(filename2,'r')
    inFilePlus2 = file(filename3,'r')
    inFileMinus2 = file(filename4,'r')
   
    senseReads=[]
    antiReads=[]
    geneNames=[]
    ratio = []
    win=[]
    win2=[]
    all_genes = {}
    count = 0
    averageGenePosPlus={}
    averageGenePosMinus={}
    averageGenePlus=[]
    averageGeneMinus=[]
    averageGenePosPlus2={}
    averageGenePosMinus2={}
    averageGenePlus2=[]
    averageGeneMinus2=[]
    geneCounts={}
    geneRejects=[]
    geneNum=0
    for i,line in enumerate(open(tssFile)):
        
        pairinfo = line.replace('\n','').split('\t')
        if i%1000 ==0:
            print i
        
        start = int(pairinfo[1])
        stop = int(pairinfo[7])
        gene1_name = pairinfo[4]
        gene1_length = int(pairinfo[2])-start
        gene2_length = stop-int(pairinfo[6])
        if gene1_length>minGenelength and gene2_length>minGenelength:
            pair_reads_plus = getGeneReads(inFilePlus,startIndex[gene1_name][0],start,stop)
            pair_reads_minus = getGeneReads(inFileMinus,startIndex[gene1_name][1],start,stop)
            pair_reads_plus2 = getGeneReads(inFilePlus2,startIndex2[gene1_name][0],start,stop)
            pair_reads_minus2 = getGeneReads(inFileMinus2,startIndex2[gene1_name][1],start,stop)
        else: continue
        if pairinfo[3]=='-':
            startp=stop-int(pairinfo[2])
            pair_reads_plus.reverse()
            pair_reads_minus.reverse()
            temp=pair_reads_plus
            pair_reads_plus=pair_reads_minus
            pair_reads_minus=temp
            pair_reads_plus2.reverse()
            pair_reads_minus2.reverse()
            temp=pair_reads_plus2
            pair_reads_plus2=pair_reads_minus2
            pair_reads_minus2=temp
        else:
            startp=int(pairinfo[6])-start
               
        norm = float(sum(pair_reads_minus[startp-600:startp]))/600
        norm2 = float(sum(pair_reads_minus2[startp-600:startp]))/600
    
        if norm>lower and norm<upper and norm2>lower and norm2<upper:
                if not normalize:norm=1
                geneNum+=1
                for index,read in enumerate(pair_reads_minus):
                    readPlus=pair_reads_plus[index]
                    pos=index-startp
                    if pos in averageGenePosPlus:
                        averageGenePosMinus[pos]+=read/norm
                        averageGenePosPlus[pos]+=readPlus/norm
                        geneCounts[pos]+=1
                    else:
                        averageGenePosMinus[pos]=read/norm
                        averageGenePosPlus[pos]=readPlus/norm
                        geneCounts[pos] = 1
                for index,read in enumerate(pair_reads_minus2):
                    readPlus=pair_reads_plus2[index]
                    pos=index-startp
                    if pos in averageGenePosPlus2:
                        averageGenePosMinus2[pos]+=read/norm
                        averageGenePosPlus2[pos]+=readPlus/norm
                        
                    else:
                        averageGenePosMinus2[pos]=read/norm
                        averageGenePosPlus2[pos]=readPlus/norm
                        
                
    positions = averageGenePosPlus.keys()
    positions.sort()
    for pos in positions:
        averageGenePlus.append(float(averageGenePosPlus[pos])/geneCounts[pos])
        averageGeneMinus.append(float(averageGenePosMinus[pos])/geneCounts[pos])
    positions2 = averageGenePosPlus.keys()
    positions2.sort()
    for pos in positions2:
        averageGenePlus2.append(float(averageGenePosPlus2[pos])/geneCounts[pos])
        averageGeneMinus2.append(float(averageGenePosMinus2[pos])/geneCounts[pos])
    
        
    print geneNum
    return positions,averageGenePlus,averageGeneMinus,positions2,averageGenePlus2,averageGeneMinus2
def runningAverage(data,window):
    genelength = len(data)
    win_starts = [0]*(window/2)
    win_starts+=range( genelength-(window) )
    win_starts+=[genelength-window]*(window/2)
    #win_stops=[window]*(window/2)
    
    win_stops = range(window/2,genelength)
    win_stops+=[genelength]*(window/2)
    norms = [np.mean(data[win_starts[index]:win_stops[index]]) for index,datum in enumerate(data)]
    return norms
def getAllGenes(inFileBase,length,lower=0,upper = 100000, normalize=1):
    inFilePlus = open(inFileBase+'_plus.txt','r')
    inFileMinus = open(inFileBase+'_minus.txt','r')
    allGenes={}
    geneLengthMax=0
    geneNum=0
    for i,line in enumerate(open(inFileBase+'_index.txt')):
        if i%1000==0:print i
        fields=line.replace('\n','').split('\t')
        strand=fields[4]
        start=int(fields[2])
        stop=int(fields[3])
        if stop-start>=length:
            if stop-start>geneLengthMax:
                geneLengthMax=stop-start
            if strand == '+':
                startIndex=int(fields[5])
                geneReads = getGeneReads(inFilePlus,startIndex,start,stop)
            else:
                startIndex=int(fields[6])
                geneReads = getGeneReads(inFileMinus,startIndex,start,stop)
                geneReads.reverse()
            norm = float(sum(geneReads[100:500]))/400
            
            if norm>lower and norm<upper:
                if not normalize:norm=1
                geneNum+=1
                allGenes[fields[0]]=[geneReads,norm]
                
            
    print geneNum
    return allGenes,geneLengthMax



