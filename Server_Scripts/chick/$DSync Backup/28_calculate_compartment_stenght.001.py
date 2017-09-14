import sys
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData

from mirnylib.numutils import observedOverExpected

from mirnylib.systemutils import setExceptionHook
setExceptionHook()

if len(sys.argv < 3):
	print "Usage: ",sys.argv[0]+" heatmap_file E1_file"

hmap = sys.argv[1]
E1_file = sys.argv[2]

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
import figPath
import ntpath
figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"_"+'Compartment_strength'


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

resolution = int(heatmap_filepath.split("-")[-1].split("k")[0])*1000
print "Resolution determined: ",resolution
print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad(heatmap_filepath, 'heatmap')
q=BD.dataDict['heatmap']
obs_exp = observedOverExpected(q)

E1_values = np.loadtxt(E1_file)
assert E1_values.shape()[0] = genome_db.numBins
bins = dict([(i,np.zeros(len(genome_db.chrmLenBins[i]))) for i in range(genome_db.chmCount)])

for i in E1_values:
	chr = genome_db.label2idx[i["f1"]]
	nt_start = genome_db.label2idx[i["f2"]]
	bin_start = nt_start/resolution
	assert bins[chr][bin_start] == 0
	bins[chr][bin_start] = i

