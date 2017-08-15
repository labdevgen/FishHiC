import sys
import os
try: 
	import numpy as np 
except:
	sys.path = ["/usr/lib64/python2.7/site-packages"]+sys.path
	import numpy as np
print "Numpy inported!"
#from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome
from mirnylib import h5dict
from getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName
setExceptionHook()

genomeName = "mm10"
genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/"+genomeName,readChrms=["#","X","Y"])


#hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res25k_hiRes.hm.IC.cis.hdf5"
hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res25k_hiRes.hm.IC.cis.hdf5"
#hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5"
resolution = extractResolutionFromFileName(hm_file)


distances = range(2,11)+[30,50]
min_distance_from_chrm_end_in_bins = 100000/resolution

result = {}
for d in distances:
	result[d] = []

for data in get_chromosomes(hm_file,genome_db,resolution):
	for dist in distances:
		for bin in xrange(dist+min_distance_from_chrm_end_in_bins,len(data)-dist-min_distance_from_chrm_end_in_bins):
			left = data[bin][bin-dist]
			right = data[bin][bin+dist]
			if left*right == 0:
				continue
			ratio = max(left,right) / min (left,right)
			result[dist].append(ratio)
print "min\tmax\tmean\tmedian\n"

if os.path.exists(hm_file+".stat"):
	import shutil
	shutil.rmtree(hm_file+".stat")

os.mkdir(hm_file+".stat")

for dist in result:
	header = "Ratio of contacts at a distance "+str(dist)+" bins for each bin at a resolution "+str(resolution)
	np.savetxt(hm_file+".stat/res"+str(resolution/1000)+"k.dist_"+str(dist)+"_bins.dist.txt",np.array(result[dist]),header=header)
	print np.min(result[dist]),np.max(result[dist]),np.mean(result[dist]),np.median(result[dist])
