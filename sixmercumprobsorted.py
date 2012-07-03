import sys
import scipy.stats
import numpy
 
def sixmercumprobsorted(sixmer_filename,pause_filename,out_filename):
  
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
       
  consensus_count=0
  halfconsensus_count=0
  noconsensus_count=0
  consensus='GGTTGC'
  for line in sixmer_prob_list:
    sixmer=line[0]
    #i=0
    #match_count=0
    #mismatch_count=0
    minus_one=line[0][-1]
    minus_two=line[0][-2]
    minus_three=line[0][-3]
    minus_four=line[0][-4]
    minus_ten=line[0][-5]
    minus_eleven=line[0][-6]
    if (minus_eleven=='G' or minus_ten=='G') and (minus_four=='T'or minus_four=='C') and minus_three=='T' and minus_two=='G' and (minus_one=='C' or minus_one=='T'):
      print 'match',sixmer
      consensus_count+=float(line[1])
    elif (minus_eleven=='G' or minus_ten=='G') and (minus_one=='C' or minus_one=='T'):
      #print 'half match',sixmer
      halfconsensus_count+=float(line[1])
    else:
      #print 'no match',sixmer
      noconsensus_count+=float(line[1])
    #print sixmer,minus_one,minus_two,minus_three,minus_four,minus_ten,minus_eleven
    #for base in sixmer:
    #  if base==consensus[i]:
    #    match_count+=1
    #  elif base!=consensus[i]:
    #    mismatch_count+=1
    #  i+=1
    #if match_count>=4:
    #  print 'match',sixmer
    #  consensus_count+=float(line[1])
    #elif match_count<4 and match_count>1:
    #  print 'half match',sixmer
    #  halfconsensus_count+=float(line[1])
    #elif match_count<2:
    #  print 'no match',sixmer
    #  noconsensus_count+=float(line[1])
      
  out_file.write(str(consensus_count)+' '+str(halfconsensus_count)+' '+str(noconsensus_count))
    
  sixmer_file.close()
  pause_file.close()
  out_file.close()
  
def main():
  sixmercumprobsorted(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == '__main__':
  main()