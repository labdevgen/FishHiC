import os
import numpy as np
import sys

from mirnylib import genome
from mirnylib import h5dict
from hiclib import binnedData
from scipy.stats import pearsonr

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
import figPath
figurePath = figPath.figure_path+"correlations_res"

from mirnylib.systemutils import setExceptionHook
setExceptionHook()


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

def binwise_pearson_correlation(q1,q2):
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
	files = [base_folder+"/"+i for i in files if i.endswith(".hm")]
	resolutions=[100,1000]
	files = [i for i in files if int(i.split("-")[-1].split("k")[0]) in resolutions]
	
print "Processing following files (",len(files),"):\n","\n".join(files)
data_list={}
for file in files:
	raw_heatmap = h5dict.h5dict(file, mode='r') 
	resolution = int(raw_heatmap['resolution'])
	BD = binnedData.binnedData(resolution, genome_db)
	BD.simpleLoad(file, 'HindIII')
	q=BD.dataDict['HindIII']
	if resolution in data_list:
		data_list[resolution][file.split("/")[-1]]=q
	else:
		data_list[resolution] = {}
		data_list[resolution][file.split("/")[-1]]=q
	print len(data_list[resolution].keys()),data_list[resolution].keys()
	
print "Calculating Pearson correlation"
results = {}

for res in sorted(data_list):
	results[res]=np.zeros(shape=(len(data_list[res]),len(data_list[res])))
	genome_db.setResolution(resolution)
	for ind_i,i in enumerate(sorted(data_list[res])[:-1]):
		for ind_j,j in enumerate(sorted(data_list[res])[ind_i+1:]):
			if i != j:
				b = binwise_pearson_correlation(data_list[res][i],data_list[res][j])
				results[res][ind_i,ind_i+1+ind_j] = b
				
for res in sorted(data_list):
	np.savetxt(figurePath+str(res)+".txt",results[res],header="\t".join(sorted(data_list[res])))