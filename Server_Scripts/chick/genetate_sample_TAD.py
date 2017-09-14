import numpy as np
import math
def fillTAD(st,Tlen,base,a,const=0):
  for i in xrange(Tlen):
    idx1 = np.array(range(0,Tlen-i))+st
    idx2 = np.array(range(i,Tlen))+st
    a[idx1,idx2] = math.pow(base,(len(a)-i))+const
    a[idx2,idx1] = math.pow(base,(len(a)-i))+const
  
a=np.zeros(shape=(100,100))
fillTAD(0,100,1.1,a)
for st in range(0,50,15):
  fillTAD(st,10,1.12,a,const=100)

fillTAD(95,5,1.13,a,const=300)
np.savetxt("chrT.txt",a)