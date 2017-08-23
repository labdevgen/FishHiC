import os
import numpy as np
import sys

from mirnylib import genome
from mirnylib import h5dict
from hiclib import binnedData
from scipy.stats import pearsonr

genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

def binwise_pearson_correlation(q1,q2,resolution):
	genome_db.setResolution(resolution)
	corr=[]
	for i in xrange(genome_db.chrmCount):
			corr+=[pearsonr(k,m)[0] for k,m in zip(q1[genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i],genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i]],
											q2[genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i],genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i]])]
	corr = np.array(corr)
	corr = corr[np.logical_not(np.isnan(corr))]
	return np.average(corr)
	
if len(sys.argv)>1:
	files=sys.argv[1:]
else:
#	base_folder='/mnt/storage/home/vsfishman/HiC/data/chick/mapped-galGal4_1MB/galGal4_1MB'
	base_folder='/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/'
	files = os.listdir(base_folder)
	files = [i for i in files if i.endswith(".hm")]
	min_resolution=40
	files = [i for i in files if int(i.split("-")[-1].split("k")[0])>=min_resolution]
	
print "Processing following files:","\n".join(files)
data_list={}
for file in files:
	file = base_folder+"/"+file
	raw_heatmap = h5dict.h5dict(file, mode='r') 
	resolution = int(raw_heatmap['resolution'])
	BD = binnedData.binnedData(resolution, genome_db)
	BD.simpleLoad(file, 'HindIII')
	q=BD.dataDict['HindIII']
	data_list[file.split("/")[-1]]=(q,resolution)
	
print "Calculating Pearson correlation"
l=len(sorted(data_list.keys()))
a=np.zeros(shape=(l,l))
for ind,val in enumerate(sorted(data_list.keys())):
	for j_ind,j_val in enumerate(sorted(data_list.keys())[ind+1:]):
		if data_list[j_val][1] != data_list[val][1]:
			print "Skipping ",val,j_val," : different resolution"
		else:
			b = binwise_pearson_correlation(data_list[val][0],data_list[j_val][0],data_list[j_val][1])
			print val,j_val,b
			a[ind,j_ind] = b
			a[j_ind,ind] = b
print a