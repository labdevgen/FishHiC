#Please note
#This code is mainly based on the doSaddle function from singleShared.py
#Which is a part of Mirny Lab HiCLib libraries, version 85979ac
import sys
import numpy as np
from scipy.stats.stats import spearmanr

import os

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

from mirnylib.numutils import EIG

if len(sys.argv) < 3:
	print "Usage: ",sys.argv[0]+" heatmap_file E1_file"
	sys.exit()

heatmap_filepath = sys.argv[1]
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

E1_values = np.genfromtxt(E1_file,dtype=None)['f2']
assert len(E1_values)==len(q)

saddles={}
for chrom in range(genome_db.chrmCount):
	saddle = np.ones((5,5), dtype = float)
	st = genome_db.chrmStartsBinCont[chrom]
	end = genome_db.chrmEndsBinCont[chrom]
	cur = q[st:end,st:end]
	
	E1 = E1_values[st:end]

	mask = np.sum(cur , axis=0) > 0
	if sum(mask)>5:
			cur = cur [mask]
			cur = cur [:,mask]
			
			cur = observedOverExpected(cur)
			E1 = E1[mask]
			assert cur.shape[0]==cur.shape[1]==len(E1)
			
			for i in range(5):
				for j in range(5):
					P1, P2 = np.percentile(E1, [20 * i, 20 * i + 20])
					mask1 = (E1 > P1) * (E1 < P2)
					P1, P2 = np.percentile(E1, [20 * j, 20 * j + 20])
					mask2 = (E1 > P1) * (E1 < P2)
					if sum(mask1)*sum(mask2) != 0:
						saddle[i, j] = np.nanmean(cur[np.ix_(mask1, mask2)])
					else:
						saddle[i, j] = None
			saddles[genome_db.idx2label[chrom]]=saddle
	else:
		print "Ommiting chromsome ",genome_db.idx2label[chrom]

vmin = min([np.min(i) for i in saddles.values()])
vmax = max([np.max(i) for i in saddles.values()])
		
if not os.path.isdir(figure_path):
	os.mkdir(figure_path)

figure_path+="/"

print "Saving results to",figure_path

for i in saddles:
	print "saving ",i
	plt.clf()
	np.savetxt(figure_path+i+".txt",saddles[i])
	plotting.plot_matrix(np.log(saddles[i]))
	plt.savefig(figure_path+i+".png",dpi=300)

all_average=np.zeros((5,5),dtype=float)
for i in range(5):
	for j in range(5):
		all_average[i,j] = np.average([saddles[c][i,j] for c in saddles if not np.isnan(saddles[c][i,j])])

all_mean=np.zeros((5,5),dtype=float)
for i in range(5):
	for j in range(5):
		all_mean[i,j] = np.nanmean([saddles[c][i,j] for c in saddles if not np.isnan(saddles[c][i,j])])

plt.clf()
np.savetxt(figure_path+"all_average.txt",all_average)
#plotting.plot_matrix(all,vmin=vmin,vmax=vmax)
plotting.plot_matrix(np.log(all_average))
plt.savefig(figure_path+"all_average.png",dpi=300)

plt.clf()
np.savetxt(figure_path+"all_mean.txt",all_mean)
#plotting.plot_matrix(all,vmin=vmin,vmax=vmax)
plotting.plot_matrix(np.log(all_mean))
plt.savefig(figure_path+"all_mean.png",dpi=300)