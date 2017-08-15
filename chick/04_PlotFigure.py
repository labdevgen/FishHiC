import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

genomeName = "galGal4_1MB"
#genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/galGal4/"+genomeName,readChrms=[],chrmFileTemplate="chr%s.fa",gapFile="galGal4.gap.txt")

genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

#base_folder='/mnt/storage/home/vsfishman/HiC/data/chick/mapped-galGal4_1MB/newData/galGal4_1MB'
base_folder='/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/'

files = os.listdir(base_folder)
files = [i for i in files if i.endswith("1000k.hm")]
min_resolution=200
files = [i for i in files if int(i.split("-")[-1].split("k")[0])>min_resolution]
print "Processing following files:","\n".join(files)
for file in files:
	file = base_folder+"/"+file
	figure_path=file+'.png'
	raw_heatmap = h5dict.h5dict(file, mode='r') 
	resolution = int(raw_heatmap['resolution'])
	BD = binnedData.binnedData(resolution, genome_db)
	BD.simpleLoad(file, 'HindIII')
	BD.removeBySequencedCount(0.5)
	# Remove 1% of regions with low coverage.
	BD.removePoorRegions(cutoff=1)
	# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
	BD.truncTrans(high=0.0005)
	# Perform iterative correction.
	BD.iterativeCorrectWithoutSS()
	
	q=np.log(BD.dataDict['HindIII'])
	#if domain_res < 1000000:
	#	print "Matrix is too big, have to resize it"
	#	q=resize_matrix(q,1000000/domain_res,np.max)
	print "saving ",figure_path
	plotting.plot_matrix(q)
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=600)
	f.close()
	plt.clf()