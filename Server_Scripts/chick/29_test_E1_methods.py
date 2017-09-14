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
from matplotlib import gridspec

heatmap_filepath = sys.argv[1]

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
import figPath
import ntpath
figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"_"+'Compartment_strength'


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

def get_by_chr_E1(genome_db,resolution):
	if heatmap_filepath.endswith(".IC"):
		raw = heatmap_filepath[:-3]
	else:
		raw = heatmap_filepath
	
	print "Using raw heatmap ",raw
	global BD_raw 
	BD_raw = binnedData.binnedData(resolution, genome_db)
	BD_raw.simpleLoad(raw, 'heatmap')
	BD_raw.removeDiagonal()

	# Remove bins with less than half of a bin sequenced
	BD_raw.removeBySequencedCount(0.5)
	# We'll do iterative correction and Eigenvector expansion on trans data only!
	# We want to remove cis, because later we want to remove poor regions in trans
	BD_raw.removeCis()
	# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts)
	# Do this before removing poor regions, because single blowouts may give
	# lots of contacts to a region which does not have much contacts otehrwise.
	BD_raw.truncTrans(high=0.0005)
	# Remove 1% of regions with low coverage
	BD_raw.removePoorRegions(cutoff=1)
	# Fake cis counts. Data gets iteratively corrected during this process...
	BD_raw.fakeCis()
	# Remove bins with zero counts for eigenvector analysis --> This will be done for each chromosome in for loop
#	BD.removeZeros()
	
	# Perform eigenvector expansion.
	

	result={"OE":{},"Classic":{},"genome_wide_Classic":{}}
	genom_wide_E1 = np.genfromtxt(raw+".eig",dtype=None)['f2']
	for chrom in range(genome_db.chrmCount):
		st = genome_db.chrmStartsBinCont[chrom]
		end = genome_db.chrmEndsBinCont[chrom]
		cur = BD_raw.dataDict['heatmap'][st:end,st:end]
		mask = np.sum(cur , axis=0) > 0
		if sum(mask) > 5:
			cur = cur [mask]
			cur = cur [:, mask]
			currentEIG, eigenvalues = EIG(cur, numPCs=1)
			if spearmanr(currentEIG[0],BD_raw.trackDict["GC"][st:end][mask])[0] < 0:
				currentEIG[0] = -currentEIG[0]
			E1=np.empty(shape=(len(mask),))*np.nan
			E1[mask] = currentEIG[0]
			result["Classic"][chrom] = E1
			
			cur = observedOverExpected(cur)
			mask = np.sum(cur , axis=0) > 0
			if sum(mask) > 5:
				cur = cur [mask]
				cur = cur [:, mask]
				currentEIG, eigenvalues = EIG(cur, numPCs=1)
				if spearmanr(currentEIG[0],BD_raw.trackDict["GC"][st:end][mask])[0] < 0:
					currentEIG[0] = -currentEIG[0]
				E1=np.empty(shape=(len(mask),))*np.nan
				E1[mask] = currentEIG[0]
				result["OE"][chrom] = E1
	
		result["genome_wide_Classic"][chrom] = genom_wide_E1[st:end]
	return result


resolution = int(heatmap_filepath.split("-")[-1].split("k")[0])*1000
print "Resolution determined: ",resolution
print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad(heatmap_filepath, 'heatmap')

q=BD.dataDict['heatmap']
E1_values = get_by_chr_E1(genome_db,resolution)
chrom = genome_db.label2idx["chr5"]
st = genome_db.chrmStartsBinCont[chrom]
end = genome_db.chrmEndsBinCont[chrom]

#fig = plt.figure(figsize=(40, 43)) 
gs = gridspec.GridSpec(4, 1, height_ratios=[15,1,1,1]) 
ax0=plt.subplot(gs[0])
ax0.imshow(np.log(q[st:end,st:end]),interpolation='nearest')
ax0.set_aspect('auto')
plt.subplot(gs[1],sharex=ax0).plot(range(len(E1_values["OE"][chrom])),E1_values["OE"][chrom])
plt.subplot(gs[2],sharex=ax0).plot(range(len(E1_values["Classic"][chrom])),E1_values["Classic"][chrom])
plt.subplot(gs[3],sharex=ax0).plot(range(len(E1_values["genome_wide_Classic"][chrom])),E1_values["genome_wide_Classic"][chrom])

plt.tight_layout()

plt.savefig("TestE1.png")