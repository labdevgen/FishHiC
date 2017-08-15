import numpy as np
import os
import gzip

import math
import subprocess
from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC
from mirnylib.numutils import fillDiagonal

domain_res='fragment'

genome_name='mm9'
genome_folder='/mnt/storage/home/vsfishman/HiC/fasta/'
genome_db = genome.Genome(genome_folder+genome_name, readChrms=['#', 'X'])

base_filename="Sp_full"
base_folder='/mnt/storage/home/vsfishman/HiC/data/'

maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'

base_out_folder = "/mnt/storage/home/vsfishman/tmp/HiC_tmp/data/"

B = 5.0
c = 1.12

def DistancetoBinN(distance):
#	return math.ceil(math.log(distance/B,c)))
	return int(distance/B)

def BinNtoDistance(n):
#	return B*math.pow(c,n)
	return (n+1)*B

def Generate_whole_genome_chromosome_file():
	max_dist_kb = 100.0
	min_dist_kb = 5.0


	fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'
	o_file = base_out_folder + "Ps_distr/"+base_filename

	if not os.path.isfile(fragment_dataset_filename):
		fragments = fragmentHiC.HiCdataset(
			filename=fragment_dataset_filename,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='w')
		fragments.parseInputData(
              dictLike=maped_reads_filepath, removeSS=True)

		fragments.filterRsiteStart(offset=5)
		fragments.filterDuplicates()
		fragments.filterLarge()
		fragments.filterExtreme(cutH=0.005, cutL=0)
	else:
		print "fragment_dataset "+ fragment_dataset_filename + " found.\n IMPORTANT: considering all requerd filters have been done"
		fragments = fragmentHiC.HiCdataset(
			filename=fragment_dataset_filename,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='a')

	Ps = np.zeros(DistancetoBinN(max_dist_kb)+1)
	
	distances = fragments.distances
	l = len(distances)
	
	count = 0
	for i in xrange(l):
		if (i % (l/10)) == 0:
			print (i / (l/10)),"0 %"
		d = distances[i]
		if d == -1:
			continue
		count += 1
		d = abs(d)/1000.0
		if (d > max_dist_kb) or ((d < min_dist_kb)):
				continue
		
		cell = DistancetoBinN(d)
		Ps[cell] += 1
			
	Ps = Ps / count

	f_out = open (o_file+".ps","w")
	for i in xrange(len(Ps)):
		f_out.write(str(BinNtoDistance(i)))
		f_out.write("\t")
		f_out.write(str(Ps[i]))
		f_out.write("\n")
	f_out.close()


base_names=["Sp_full","Fib_full","ESC_full"]
#base_names=["Sp_1000","ESC_full","Fib_full"]
#base_names=["Sp_1000"]
for j in base_names:
		base_filename=j
		Generate_whole_genome_chromosome_file()