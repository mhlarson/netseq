import sys
import scipy.stats
import numpy
 
def sixmercumprob(sixmer_filename,pause_filename,out_filename):
  
  sixmer_file = open(sixmer_filename, 'rU')
  sixmer_string = sixmer_file.readlines()
  sixmer_list=[]
  sixmer_prob_list=[]
  
  pause_file = open(pause_filename, 'rU')  ## clean pause file
  pause_string = pause_file.readlines()
  pause_dict={}
  
  out_file = open(out_filename,'wb')  
  
  for line in sixmer_string:
    sixmer=line.strip('\n').split('\t')[0]
    sixmer_list.append(sixmer)
      
  for line in pause_string:
    pause=line.strip('\n')
    pause_sixmer=pause[4]+pause[5]+pause[11]+pause[12]+pause[13]+pause[14]
    if pause_sixmer in pause_dict:
      pause_dict[pause_sixmer]+=1
    else:
      pause_dict[pause_sixmer]=1
      
  for line in sixmer_list:
    sixmer=line
    if sixmer in pause_dict:
      sixmer_prob_list.append((sixmer,pause_dict[sixmer]))
    else:
      sixmer_prob_list.append((sixmer,0))
       
  count=0  
  for line in sixmer_prob_list:
    sixmer=line[0]
    count+=float(line[1])
    out_file.write(sixmer+'\t'+str(count)+'\n')
    
  sixmer_file.close()
  pause_file.close()
  out_file.close()
  
def main():
  sixmercumprob(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()