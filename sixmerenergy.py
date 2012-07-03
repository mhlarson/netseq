import sys
import scipy.stats
import numpy
 
def sixmerenergy(sixmer_filename,directoryname,out_filename):
  
  sixmer_file = open(sixmer_filename, 'rU')
  sixmer_string = sixmer_file.readlines()
  sixmer_dict={}
  
  out_file = open(out_filename,'wb')  
  
  for line in sixmer_string:
    sixmer=line.strip('\n')
    sixmer_dict[sixmer]=0
  
  A_file = open(directoryname+'G_A.csv','rU')
  C_file = open(directoryname+'G_C.csv','rU')
  G_file = open(directoryname+'G_G.csv','rU')
  T_file = open(directoryname+'G_T.csv','rU')
  
  A_string = A_file.readlines()
  A_dict={}
  j=19
  for line in A_string:
    A_energy=float(line.strip('\n'))
    A_dict[j]=A_energy
    j+=-1
    
  C_string = C_file.readlines()
  C_dict={}
  j=19
  for line in C_string:
    C_energy=float(line.strip('\n'))
    C_dict[j]=C_energy
    j+=-1
      
  G_string = G_file.readlines()
  G_dict={}
  j=19
  for line in G_string:
    G_energy=float(line.strip('\n'))
    G_dict[j]=G_energy
    j+=-1
  
  T_string = T_file.readlines()
  T_dict={}
  j=19
  for line in T_string:
    T_energy=float(line.strip('\n'))
    T_dict[j]=T_energy
    j+=-1
      
  for line in sixmer_dict:
    sixmer=line
    n_11=line[0]
    n_10=line[1]
    n_4=line[2]
    n_3=line[3]
    n_2=line[4]
    n_1=line[5]
    if n_11=='A':
      n_11_energy=A_dict[-11]
    if n_11=='C':
      n_11_energy=C_dict[-11]
    if n_11=='G':
      n_11_energy=G_dict[-11]
    if n_11=='T':
      n_11_energy=T_dict[-11]
    
    if n_10=='A':
      n_10_energy=A_dict[-10]
    if n_10=='C':
      n_10_energy=C_dict[-10]
    if n_10=='G':
      n_10_energy=G_dict[-10]
    if n_10=='T':
      n_10_energy=T_dict[-10]
    
    if n_4=='A':
      n_4_energy=A_dict[-4]
    if n_4=='C':
      n_4_energy=C_dict[-4]
    if n_4=='G':
      n_4_energy=G_dict[-4]
    if n_4=='T':
      n_4_energy=T_dict[-4]
      
    if n_3=='A':
      n_3_energy=A_dict[-3]
    if n_3=='C':
      n_3_energy=C_dict[-3]
    if n_3=='G':
      n_3_energy=G_dict[-3]
    if n_3=='T':
      n_3_energy=T_dict[-3]
    
    if n_2=='A':
      n_2_energy=A_dict[-2]
    if n_2=='C':
      n_2_energy=C_dict[-2]
    if n_2=='G':
      n_2_energy=G_dict[-2]
    if n_2=='T':
      n_2_energy=T_dict[-2]
    
    if n_1=='A':
      n_1_energy=A_dict[-1]
    if n_1=='C':
      n_1_energy=C_dict[-1]
    if n_1=='G':
      n_1_energy=G_dict[-1]
    if n_1=='T':
      n_1_energy=T_dict[-1]
    
    sixmer_energy = n_11_energy + n_10_energy + n_4_energy + n_3_energy + n_2_energy + n_1_energy
    sixmer_dict[sixmer]=sixmer_energy
    
  for line in sixmer_dict:
    sixmer=line
    out_file.write(sixmer+'\t'+str(sixmer_dict[sixmer])+'\n')
  
  sixmer_file.close()
  out_file.close()
  A_file.close()
  C_file.close()
  G_file.close()
  T_file.close()
  
def main():
  sixmerenergy(sys.argv[1],sys.argv[2],sys.argv[3])
  
if __name__ == '__main__':
  main()