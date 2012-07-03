import sys
import scipy.stats
import numpy
 
def PearsonCorrelation(wave1_filename,wave2_filename):
  
  wave1_file = open(wave1_filename, 'rU')
  wave2_file = open(wave2_filename, 'rU')
  
  wave1_string = wave1_file.readlines()
  wave2_string = wave2_file.readlines()
  
  wave1_list=[]
  wave2_list=[] 
  variance1_list=[]
  variance2_list=[]
  std1_list=[]
  std2_list=[]
  weightedwave1_list=[]
  weightedwave2_list=[]
  
  covar_list=[]
  
  weight1=1#/(0.02)
  weight2=1#/(0.08)
  
  for line in wave1_string:
    value1=float(line.strip('\n'))
    wave1_list.append(value1)
    weightedvalue1=value1*weight1
    weightedwave1_list.append(weightedvalue1)
  
  for line in wave2_string:
    value2=float(line.strip('\n'))
    wave2_list.append(value2)
    weightedvalue2=value2*weight2
    weightedwave2_list.append(weightedvalue2)
  
  weightedmean1=numpy.sum(weightedwave1_list)/(1*weight1)
  weightedmean2=numpy.sum(weightedwave2_list)/(1*weight2)
  
  mean1=numpy.mean(wave1_list)
  mean2=numpy.mean(wave2_list)
  
  for element1 in wave1_list:
    #variance1=element1-mean1
    variance1=weight1*(element1-weightedmean1)
    variance1_list.append(float(variance1))
    std1=(variance1)**2
    std1_list.append(float(std1))
   
  for element2 in wave2_list:
    #variance2=element2-mean2
    variance2=weight2*(element2-weightedmean2)
    variance2_list.append(variance2)
    std2=(variance2)**2
    std2_list.append(float(std2))
  
  for i in range(0,160):
    covar=weight1*(variance1_list[i] * variance2_list[i])
    covar_list.append(covar)
    
  covarsum=(numpy.sum(covar_list))/(1*weight1)
  
  sumstd1=(numpy.sum(std1_list))**0.5
  sumstd2=(numpy.sum(std2_list))**0.5
  
  r=covarsum/(sumstd1*sumstd2)
  #print r
  
  r=scipy.stats.pearsonr(wave1_list,wave2_list)
  print r
  
  wave1_file.close()
  wave2_file.close()
  
def main():
  PearsonCorrelation(sys.argv[1],sys.argv[2])

if __name__ == '__main__':
  main()