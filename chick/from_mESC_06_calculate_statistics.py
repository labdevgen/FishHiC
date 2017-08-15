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

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")
from getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName

setExceptionHook()


genome_db_contig = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")


hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII-100k.hm.IC"
resolution = extractResolutionFromFileName(hm_file)

distances = range(10,50,2)
min_distance_from_chrm_end_in_bins = 100000/resolution

result = {}
for d in distances:
	result[d] = []

print "Calculating..."
for data in get_chromosomes(hm_file,genome_db_contig,resolution):
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
print "Done!"