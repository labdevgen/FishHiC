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
from mirnylib import h5dict

setExceptionHook()

import os
#from defineGenome import getGenome

genomeName = "mm10"
genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/"+genomeName,readChrms=["#","X","Y"])

data_folder = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/"
from getIntraChrHeatmaps import extractResolutionFromFileName

def filter_hires_heatmap(mode="cis",hm_file=""):
	from hiclib import highResBinnedData

	resolution = extractResolutionFromFileName(hm_file)
	if resolution == None:
		raise

	# Create a  object, load the data.
	print "creating an object"
	hmap = highResBinnedData.HiResHiC(genome_db,resolution)
			
	print "loading data"
	hmap.loadData(hm_file, mode=mode)

	print "saving pict of heatmap"
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from mirnylib import plotting
	
	chr0array = hmap.data[(0,0)].getData()
	maxlen=min(10000,len(chr0array))

	to_plot=chr0array[0:maxlen,0:maxlen]
	figure_path = hm_file+"stage1.png"
	print "saving ",figure_path
	plotting.plot_matrix(np.log(to_plot))
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=300)
	f.close()
	plt.clf()

	
	# Remove the contacts between loci located within the same bin +/- 1 bin.
	hmap.removeDiagonal(m=1)

	to_plot=hmap.data[(0,0)].getData()[0:maxlen,0:maxlen]
	figure_path = hm_file+"stage2.png"
	print "saving ",figure_path
	plotting.plot_matrix(np.log(to_plot))
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=300)
	f.close()
	plt.clf()

	# Removes 0.5 percent of regions with low coverage.
	hmap.removePoorRegions(percent=0.5)

	# Perform iterative correction.
	hmap.iterativeCorrection()

	to_plot=hmap.data[(0,0)].getData()[0:maxlen,0:maxlen]
	figure_path = hm_file+"stage3.png"
	print "saving ",figure_path
	plotting.plot_matrix(np.log(to_plot))
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=300)
	f.close()
	plt.clf()

	# Save the iteratively corrected heatmap.
	hmap.export(hm_file+".IC."+mode+".hdf5")

def filter_bychr_heatmap(hm_file):
	
	resolution = extractResolutionFromFileName(hm_file)
	if resolution == None:
		raise
	from hiclib import binnedData 
	# Create a  object, load the data.
	print "creating an object"
	hmap = binnedData.binnedData(resolution,genome_db)
	
	print "loading data"
	hmap.simpleLoad(hm_file,"heatmap")
	
	print "saving pict of heatmap"
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	from mirnylib import plotting
	
	maxlen=min(10000,len(hmap.dataDict["heatmap"]))

	a=hmap.dataDict["heatmap"][0:maxlen,0:maxlen]
	figure_path = hm_file+"stage1.png"
	print "saving ",figure_path
	plotting.plot_matrix(np.log(a))
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=600)
	f.close()
	plt.clf()


	# Remove the contacts between loci located within the same bin +/- 1 bin.
	hmap.removeDiagonal(m=1)

	hmap.removeBySequencedCount()  # new filter: omit all bins with less than 0.5 coverage by sequenced bases (i.e. bases present in the genome)

	hmap.removePoorRegions(cutoff = 0.5, coverage=True)  # remove .5% bins with the lowest number of records (i.e. non-zero entrees in the matrix)
	# This filter was updated to remove bins which have zero contacts and one PCR blowout. Those bins would have many reads, but all reads will be with one or few other bins. 

	hmap.truncTrans() # remove PCR blowouts from trans data
	
	a=hmap.dataDict["heatmap"][0:maxlen,0:maxlen]
	figure_path = hm_file+"stage2.png"
	print "saving ",figure_path
	plotting.plot_matrix(np.log(a))
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=200)
	f.close()
	plt.clf()
	
	hmap.iterativeCorrectWithoutSS(force=True) #do iterative correction 
	
	a=hmap.dataDict["heatmap"][0:maxlen,0:maxlen]
	figure_path = hm_file+"stage3.png"
	print "saving ",figure_path
	plotting.plot_matrix(np.log(a))
	plt.subplots_adjust(bottom=0.15)
	f = open(figure_path, "wb")
	plt.savefig(figure_path,dpi=600)
	f.close()
	plt.clf()

	# Save the iteratively corrected heatmap.
	hmap.export("heatmap",hm_file+".IC.hdf5",False)

hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-NcoI_refined.frag_res10k_hiRes.hm"
filter_hires_heatmap(mode="cis",hm_file=hm_file)
#filter_bychr_heatmap()
