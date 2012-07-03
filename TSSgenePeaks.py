HELP_STRING = """
TSSgenePeaks.py
Takes in a track file and determines where spikes lie. Peaks are called for every
gene where the mean reads per bp are above the cutoff value. 
For every position along a gene, the mean and standard deviation is determined
for a window of size <window> surrounding that position.
The position is a pause if its value is above <peakmin> and is <fold> standard
deviations above the mean. Also the mean for the window has to be above tnorm_thresh.

The process is repeated <rounds> times. After each iteration, the peaks are
removed from the gene and ignored for mean and standard deviation calculations.
New peaks are then found.


Author: Stirling Churchman
Date: Decemeber 29, 2009
Updated: May 4, 2011

* all options overided in script! Edit script to change these.
     -h     print this help message
     -f     data file
     -s     peak strength min
     -S     peak strength max
     -c     gene levels cutoff
     -o     output file (required)
"""
 
import sys

from getopt import getopt
import analysis as a
import numpy as np
import os
import math


def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt(argv[1:], "hf:s:S:c:o:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
   
       
    for (opt, opt_arg) in optlist:
        #print opt
        #print opt_arg
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-f':
            condition = opt_arg
        elif opt == '-s':
            peakStrength = int(opt_arg)
        elif opt == '-S':
            peakStrengthMax = int(opt_arg)
        elif opt == '-c':
            cutoff= int(opt_arg)
        elif opt == '-o':
            outputFile = opt_arg
    
         
    
    condition = '/home/mattlarson/11232011/greA/greA_all'
    inputFile3 = '/home/mattlarson/genomes/MG1655_TSSs.mochiview.txt'
    genomeFile = '/home/mattlarson/genomes/NC_000913_pause.fna'
    outputBase = '/home/mattlarson/11232011/greA/pausing/'
    promoterfile = '/home/mattlarson/genomes/RegDB_allPromoters_050611.mochiview.txt'
    
    #Stirling's original settings: window=200,fold=4,peakmin=3,tnorm_thresh=2,rounds=10, peakStrengthMin=0,cutoff=0
    seqLength=15    #Length of the retrieved sequence upstream and downstream of the pause site (i.e. total sequence will be double this value)
    cutoff=1        #Genes must have a mean value larger than <cutoff> to be counted in pause analysis (i.e. threshold for total gene)
    peakStrengthMin=10 #Raw peak height divided by the mean of <window> surrounding that position (with <window>/2 on each side of the pause peak)
    window=200      #For every position, the mean and standard deviation is determined for a window of size <window> surrounding that position (i.e. <window>/2 on either side of the pause peak)
    fold=6         #Position is a pause if its value is <fold> standard deviations above the mean
    peakmin=10      #Position is a pause if its value is above <peakmin>    
    tnorm_thresh=2  #Mean for the window has to be above <tnorm_thresh> to be counted in pause analysis (i.e. treshold for the window)
    rounds=4       #Process is repeated <rounds> times. After each iteration, peaks are removed and ignored for mean and standard deviation calculations. New peaks are then found.
    dir = 'cutoff_%s_peakStrengthMin_%s_means_wo_zeros0311' %(cutoff,peakStrengthMin)
    
    if dir not in os.listdir(outputBase):
        os.mkdir(outputBase+dir)
    print condition
    print window,fold,peakmin,tnorm_thresh,rounds
    
    peakSeq,peakPos =a.getPeakSeqs(genomeFile,condition+'.plus.wig',condition+'.minus.wig',inputFile3,seqLength,cutoff,peakStrengthMin,window,fold,peakmin,tnorm_thresh,rounds,promoterfile)
    filenames=condition.split('/')
    oFile = outputBase+dir+'/'+filenames[-4]+filenames[-1]+'_%s_%s_%s_%s_%s.csv' % (window,fold,peakmin,tnorm_thresh,rounds)
    outFile = open(oFile, 'w')
    
    upstream=[seq[0] for seq in peakSeq]
    downstream=[seq[1] for seq in peakSeq]
        
    pauselogo=[]
    for uppauseindex in range(-1,-16,-1):
      A_counts=0
      C_counts=0
      G_counts=0
      T_counts=0
      for i in enumerate(upstream):
        base=i[1][uppauseindex]
        if base=='A':
          A_counts=A_counts+1
        elif base=='C':
          C_counts=C_counts+1
        elif base=='G':
          G_counts=G_counts+1
        elif base=='T':
          T_counts=T_counts+1
        
      TotalReads = (A_counts)+(C_counts)+(G_counts)+(T_counts)
      percentA=math.floor((float(A_counts)/TotalReads)*100)
      percentC=math.floor((float(C_counts)/TotalReads)*100)
      percentG=math.floor((float(G_counts)/TotalReads)*100)
      percentT=math.floor((float(T_counts)/TotalReads)*100)
      
      percentlist=['up',uppauseindex,'%A:',percentA,'%C',percentC,'%G',percentG,'%T',percentT]
      pauselogo.insert(0,percentlist)
    
    for downpauseindex in range(0,15,1):
      A_counts=0
      C_counts=0
      G_counts=0
      T_counts=0
      for i in enumerate(downstream):
        base=i[1][downpauseindex]
        if base=='A':
          A_counts=A_counts+1
        elif base=='C':
          C_counts=C_counts+1
        elif base=='G':
          G_counts=G_counts+1
        elif base=='T':
          T_counts=T_counts+1
        
      TotalReads = (A_counts)+(C_counts)+(G_counts)+(T_counts)
      percentA=math.floor((float(A_counts)/TotalReads)*100)
      percentC=math.floor((float(C_counts)/TotalReads)*100)
      percentG=math.floor((float(G_counts)/TotalReads)*100)
      percentT=math.floor((float(T_counts)/TotalReads)*100)
    
      percentlist=['down',downpauseindex+1,'%A:',percentA,'%C',percentC,'%G',percentG,'%T',percentT]
      pauselogo.append(percentlist)
    #print pauselogo
    
    #percentA = float(len([seq for seq in upstream if seq[-1]=='A']))/len(upstream)*100
    #percentT = float(len([seq for seq in downstream if seq[0]=='T']))/len(upstream)*100
    #percentATC = float(len([seq for (i,seq) in enumerate(upstream) if (seq[-1]=='A' and downstream[i][:2]=='TC') ]))/len(upstream)*100
    #print "Percent A is %s" % percentA
    #print "Percent T is %s" % percentT
    #print "Percent ATC is %s" % percentATC
    
    print "Number of pauses analyzed: %s" % len(upstream)
    print pauselogo
    
    for i,useq in enumerate(upstream):
        
        full_seq = useq+downstream[i]
        #outfile format: gene, strand, txn start, dist from txn start, pause site, peak height, peak strength, avg of reads surrounding peak, trans start, dist from trans start
        
        if peakPos[i][3]=='+':
          gene_name=peakPos[i][12]
          strand=peakPos[i][3]
          txn_start=peakPos[i][-3]
          dist_txnstart=peakPos[i][-5]
          pause_site=peakPos[i][-3]+peakPos[i][-5]
          peak_height=peakPos[i][-4]
          peak_strength=peakPos[i][-2]
          avg_gene=peakPos[i][-1]
          trans_start=peakPos[i][1]
          dist_transstart=pause_site-int(trans_start)
          fold_above_mean=math.floor(peak_height/avg_gene)
          
          outFile.write(">%s_%s_%s_%s_%s_%s_%s_%s_%s\n" % (gene_name,strand,txn_start,dist_txnstart,pause_site,peak_height,peak_strength,avg_gene,fold_above_mean))
          outFile.write("%s\n" % full_seq)
          
        elif peakPos[i][3]=='-':
          stop_tolerance=0      #Allows user to look for pauses beyond the stop boundary of the gene coding sequence
          genelength=int(peakPos[i][-3])-int(peakPos[i][1])+stop_tolerance
          gene_name=peakPos[i][12]
          strand=peakPos[i][3]
          txn_start=peakPos[i][-3]
          dist_txnstart=genelength-peakPos[i][-5]
          pause_site=txn_start-dist_txnstart
          peak_height=peakPos[i][-4]
          peak_strength=peakPos[i][-2]
          avg_gene=peakPos[i][-1]
          trans_start=peakPos[i][2]
          dist_transstart=int(trans_start)-pause_site
          fold_above_mean=math.floor(peak_height/avg_gene)
          
          outFile.write(">%s_%s_%s_%s_%s_%s_%s_%s_%s\n" % (gene_name,strand,txn_start,dist_txnstart,pause_site,peak_height,peak_strength,avg_gene,fold_above_mean))
          outFile.write("%s\n" % full_seq)
        
    
    outFile.close()
    
    
    


##############################################
if __name__ == "__main__":
    sys.exit(main())
