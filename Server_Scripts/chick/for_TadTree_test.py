import sys
import numpy as np

def fill_diag(a):
	for i in xrange(len(a)-2):
#		a[i,i] = a[i,i+2]*2
		a[i,i+1] = a[i,i+2]*1.5
		a[i+1,i] = a[i,i+2]*1.5

file = sys.argv[1]
a = np.loadtxt(file)
fill_diag(a)
#nonzero = np.nonzero(a)
#m = np.min(a[nonzero])
#a = a*(1./m)
#a = a.astype(np.int64)
np.savetxt("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/TADtree/test2.gz",a,fmt="%i",delimiter="\t")