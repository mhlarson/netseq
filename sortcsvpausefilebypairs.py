import sys
import csv

# Need to use the full pause file (i.e. with two lines for each pause)

def sortcsvpausefilebypairs(csvpausefilename):
  csvpausefile=open(csvpausefilename, 'rU')
  
  outfilename=csvpausefilename[:-4]

  outfile_GA=open(outfilename+"_GA.csv",'wb')
  outfile_GC=open(outfilename+"_GC.csv",'wb')
  outfile_GG=open(outfilename+"_GG.csv",'wb')
  outfile_GT=open(outfilename+"_GT.csv",'wb')
  
  outfile_AG=open(outfilename+"_AG.csv",'wb')
  outfile_CG=open(outfilename+"_CG.csv",'wb')
  outfile_TG=open(outfilename+"_TG.csv",'wb')
  
  #outfile_AA=open(outfilename+"_AA.csv",'wb')
  #outfile_AC=open(outfilename+"_AC.csv",'wb')
  #outfile_AG=open(outfilename+"_AG.csv",'wb')
  #outfile_AT=open(outfilename+"_AT.csv",'wb')
  
  #outfile_CA=open(outfilename+"_CA.csv",'wb')
  #outfile_CC=open(outfilename+"_CC.csv",'wb')
  #outfile_CG=open(outfilename+"_CG.csv",'wb')
  #outfile_CT=open(outfilename+"_CT.csv",'wb')
  
  #outfile_GA=open(outfilename+"_GA.csv",'wb')
  #outfile_GC=open(outfilename+"_GC.csv",'wb')
  #outfile_GG=open(outfilename+"_GG.csv",'wb')
  #outfile_GT=open(outfilename+"_GT.csv",'wb')
  
  #outfile_TA=open(outfilename+"_TA.csv",'wb')
  #outfile_TC=open(outfilename+"_TC.csv",'wb')
  #outfile_TG=open(outfilename+"_TG.csv",'wb')
  #outfile_TT=open(outfilename+"_TT.csv",'wb')
  
  #csvpausestring=csvpausefile.readlines()
  
  for line in csvpausefile:
    pause_row = line.strip('\n').split('\t')
    pause_index=pause_row[0]
    pause_sequence=pause_row[1]
    
    if pause_sequence[14:16]=='CG' or pause_sequence[14:16]=='TG':
      if pause_sequence[4:6]=='GA':
        outfile_GA.write(pause_sequence+'\n')
      elif pause_sequence[4:6]=='GC':
        outfile_GC.write(pause_sequence+'\n')
      elif pause_sequence[4:6]=='GG':
        outfile_GG.write(pause_sequence+'\n')
      elif pause_sequence[4:6]=='GT':
        outfile_GT.write(pause_sequence+'\n')
      
      elif pause_sequence[4:6]=='AG':
        outfile_AG.write(pause_sequence+'\n')
      elif pause_sequence[4:6]=='CG':
        outfile_CG.write(pause_sequence+'\n')
      elif pause_sequence[4:6]=='TG':
        outfile_TG.write(pause_sequence+'\n')
      
    #if csvpausestring[i+1][14:16]=='AA':
    #  outfile_AA.write(csvpausestring[i])
    #  outfile_AA.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='AC':
    #  outfile_AC.write(csvpausestring[i])
    #  outfile_AC.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='AG':
    #  outfile_AG.write(csvpausestring[i])
    #  outfile_AG.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='AT':
    #  outfile_AT.write(csvpausestring[i])
    #  outfile_AT.write(csvpausestring[i+1]) 
    
    #elif csvpausestring[i+1][14:16]=='CA':
    #  outfile_CA.write(csvpausestring[i])
    #  outfile_CA.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='CC':
    #  outfile_CC.write(csvpausestring[i])
    #  outfile_CC.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='CG':
    #  outfile_CG.write(csvpausestring[i])
    #  outfile_CG.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='CT':
    #  outfile_CT.write(csvpausestring[i])
    #  outfile_CT.write(csvpausestring[i+1]) 
    
    #elif csvpausestring[i+1][14:16]=='GA':
    #  outfile_GA.write(csvpausestring[i])
    #  outfile_GA.write(csvpausestring[i+1])
    #elif csvpausestring[i+1][14:16]=='GC':
    #  outfile_GC.write(csvpausestring[i])
    #  outfile_GC.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='GG':
    #  outfile_GG.write(csvpausestring[i])
    #  outfile_GG.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='GT':
    #  outfile_GT.write(csvpausestring[i])
    #  outfile_GT.write(csvpausestring[i+1]) 
    
    #elif csvpausestring[i+1][14:16]=='TA':
    #  outfile_TA.write(csvpausestring[i])
    #  outfile_TA.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='TC':
    #  outfile_TC.write(csvpausestring[i])
    #  outfile_TC.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='TG':
    #  outfile_TG.write(csvpausestring[i])
    #  outfile_TG.write(csvpausestring[i+1]) 
    #elif csvpausestring[i+1][14:16]=='TT':
    #  outfile_TT.write(csvpausestring[i])
    #  outfile_TT.write(csvpausestring[i+1]) 
    
  #outfile_AA.close()
  #outfile_AC.close()
  #outfile_AG.close()
  #outfile_AT.close()
  
  #outfile_CA.close()
  #outfile_CC.close()
  #outfile_CG.close()
  #outfile_CT.close()
  
  #outfile_GA.close()
  #outfile_GC.close()
  #outfile_GG.close()
  #outfile_GT.close()
  
  #outfile_TA.close()
  #outfile_TC.close()
  #outfile_TG.close()
  #outfile_TT.close()
  
  outfile_GA.close()
  outfile_GC.close()
  outfile_GG.close()
  outfile_GT.close()
  
  outfile_AG.close()
  outfile_CG.close()
  outfile_TG.close()
  
  csvpausefile.close()
  
def main():
  sortcsvpausefilebypairs(sys.argv[1])

if __name__ == '__main__':
  main()