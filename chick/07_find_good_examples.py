import sys
try: 
	import numpy as np 
except:
	sys.path = ["/usr/lib64/python2.7/site-packages"]+sys.path
	import numpy as np
print "Numpy inported!"
#from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome
from hiclib import highResBinnedData
from mirnylib import h5dict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
setExceptionHook()
import os

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")
from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName

from lib_coordinates_converter import *
converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
converter.prepare_acc2chrmDict()							
converter.create_agp_dict()							


genome_db_contig = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")


hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII-100k.hm.IC"

resolution = extractResolutionFromFileName(hm_file)

percentiles = [60,70,75,90,95,99,99.5,99.9,99.99,99.999]

prefered_distances= [44,46,48]
#prefered_distances= []

#find avaliable distances
if not os.path.exists(hm_file+".stat"):
	print "Distances dir %s not found" % (hm_file+".stat")
	raise "please run 06_calculate_statistics.py first"
	
files=[f for f in os.listdir(hm_file+".stat") if f.endswith(".dist.txt")]

#load statistics
distances={}
for file in files:
	res=extractResolutionFromFileName(file)
	if res != resolution:
		print "Resolution %d <> resolution of distance files %d" % (res,resolution)
		raise "Resolution mismatch"
	key = int(file.split("k.dist_")[-1].split("_")[0])
	if len(prefered_distances)==0 or (key in prefered_distances):
		distances[key]=np.loadtxt(hm_file+".stat/"+file)
	else:
		print "Skipping distance",key
print "Following distances found: ",sorted(distances.keys())
max_distance = max(distances.keys())
min_distance_from_chrm_end_in_bins = 100000/resolution


percentiles.sort()
criteias = np.zeros(shape=(len(percentiles),len(distances)))
for ind,p in enumerate(percentiles):
	for jnd,dist in enumerate(sorted(distances.keys())):
		criteias[ind,jnd] = np.percentile(distances[dist],p)


results_file_name = hm_file+".good_examples."+".".join(map(str,sorted(distances.keys())))+".txt"

np.savetxt(results_file_name+".criteias.txt",criteias,header="distances(columns): "+"\t".join(map(str,sorted(distances.keys())))+" percentiles(rows): "+ "\t".join(map(str,percentiles)))


results_file = open(results_file_name,"w")
results_file.write("Contig\tNt\tChr\tNt\tMax_percentile\t"+"\t".join(map(str,sorted(distances.keys())))+"\n")

print "Calculating..."

result_counts=dict([(p,0) for p in percentiles])

for chr_numb,data in enumerate(get_chromosomes(hm_file,genome_db_contig,resolution)):

	for bin in xrange(max_distance+min_distance_from_chrm_end_in_bins,len(data)-max_distance-min_distance_from_chrm_end_in_bins):
		max_bin_percentile = 101 #percentile can not be > 100, so 101 indicates we are in the begining of the loop
		sign = None #sign means left contacts > right contacts (True),  left contacts < right contacts (False), or not defined (None)
		ratios=[]
		for jnd,dist in enumerate(sorted(distances.keys())):
			max_curren_percentile = 101
			left = data[bin][bin-dist] #contacts on the left side of the bin
			right = data[bin][bin+dist] #contacts on the right side of the bin
			if left*right == 0: #bin has 0 contacts on left or right side
				max_bin_percentile = 101
				break #than we do not consider it at all (max_percentile will be -1)
			dist_sign = (left - right) > 0
			if sign == None:
				sign = dist_sign
			elif sign != dist_sign:
				max_bin_percentile = 101
				break			
				
			ratio=max(left,right) / min (left,right) #reatio between left and right contacts
			if left > right:
				ratios.append(ratio)
			else:
				ratios.append(-ratio)
			
			for ind,p in enumerate(percentiles): #lets check whether this ratio is higher that q-th percentile for defined distance
				if ratio > criteias[ind,jnd]:
					max_curren_percentile = p
				else: #we started from the minimal percentile, increasing the strength of criteria
					break #if we can not satisfy criteria at any point, it makes no scence to go further (break percentiles loop)
					
			if max_curren_percentile==101: #if we do not satisfy any percentile for this distance 
				max_bin_percentile = 101
				break #we just scip this bin (break distance loop)
			else:
				max_bin_percentile = min(max_bin_percentile,max_curren_percentile)
				
		if max_bin_percentile != 101: #if we satisfy any percentile for all distances
			#we need to write a result
			assert len(ratios)==len(distances)
			chr_name_contigLevel = "N"+genome_db_contig.idx2label[chr_numb]
			nt_contigLevel = int((bin*resolution)+resolution/2.)
			position_chrmLevel = converter.contig2chrm(chr_name_contigLevel,nt_contigLevel)
			results_file.write("\t".join(map(str,[chr_name_contigLevel,
													nt_contigLevel,
													position_chrmLevel.chrm,
													position_chrmLevel.nt,
													max_bin_percentile]+ratios))+"\n")
			result_counts[max_bin_percentile] += 1

results_file.close()

print "Number of bins satisfying criteias (for each percentile) are:"
for i in sorted(result_counts):
	print i,result_counts[i]
print "Results saved in \n",results_file_name