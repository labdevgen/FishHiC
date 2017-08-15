import sys
import numpy as np

eig = np.genfromtxt("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-100k.hm.eig",dtype=None)

def sortfunct(s):
		if s in ["chrW","chrZ"]:
			return 10000000
		else:
			return 100*len(s)

chromosomes = np.unique(eig["f0"])


sortorder = sorted(range(len(chromosomes)), key=lambda k: (sortfunct(chromosomes[k]),chromosomes[k]))
chromosomes=chromosomes[sortorder]


for i in chromosomes:
	e1_chr = [eig[ind]["f2"] for ind in xrange(len(eig)) if eig[ind]["f0"]==i]
	print i,np.average(e1_chr),np.median(e1_chr)

print eig
print np.average(eig["f2"]),np.std(eig["f2"]),np.median(eig["f2"])