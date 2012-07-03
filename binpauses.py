import sys
import csv

def binpauses(csvpausefilename,beforefilename,afterfilename):
  csvpausefile=open(csvpausefilename, 'rU')
    
  pauses_before_trans=[]
  pauses_after_trans=[]
  for line in csvpausefile:
    if line[0]=='>':
      dist_trans=int(line.split('_')[-1])
      pause_amplitude=[float(line.split('_')[5])]
      if dist_trans<0:
        pauses_before_trans.append(pause_amplitude)
      else:
        pauses_after_trans.append(pause_amplitude)
      
  outfilebefore=csv.writer(file(beforefilename,'wb'),dialect='excel')
  outfilebefore.writerows(pauses_before_trans)
  
  outfileafter=csv.writer(file(afterfilename,'wb'),dialect='excel')
  outfileafter.writerows(pauses_after_trans)
  
  csvpausefile.close()
  
def main():
  binpauses(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()