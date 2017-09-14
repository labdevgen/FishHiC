"""
An example script to perform iterative correction of a raw heatmap and saving eigenvectors.
"""

import os
from mirnylib import genome
from mirnylib.h5dict import h5dict
from hiclib.binnedData import binnedData
import numpy as np
import sys
sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")
from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName


genome_db_chrmLevel = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")
				
				
hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChIE-all-HindIII-100k.hm"

f_out_path=hm_file+'.eig'

NumEigenvectors = 1  # number of eigenvectors to compute

# Read resolution from one of the datasets
resolution = extractResolutionFromFileName(hm_file)

# Define the binnedData object, load data
BD = binnedData(resolution, genome_db_chrmLevel)
BD.simpleLoad(hm_file, 'heatmap')

BD.removeDiagonal()

# Remove bins with less than half of a bin sequenced
BD.removeBySequencedCount(0.5)

# We'll do iterative correction and Eigenvector expansion on trans data only!
# We want to remove cis, because later we want to remove poor regions in trans
BD.removeCis()

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts)
# Do this before removing poor regions, because single blowouts may give
# lots of contacts to a region which does not have much contacts otehrwise.
BD.truncTrans(high=0.0005)

# Remove 1% of regions with low coverage
BD.removePoorRegions(cutoff=1)

# Fake cis counts. Data gets iteratively corrected during this process...
BD.fakeCis()

# Remove bins with zero counts for eigenvector analysis
BD.removeZeros()

# Perform eigenvector expansion. Eigenvectors are now in BD.EigDict
# Corresponding eigenvalues are in BD.eigEigenvalueDict
print "Calculating Eigenvectors"
BD.doEig(numPCs=NumEigenvectors)

# BD.EigDict["HindIII"][0] is a first eigenvector of HindIII dataset
# BD.EigDict["NcoI"][2] is the third eigenvector of HindIII dataset
# etc.

# Now we restore regions with zero counts.
# We specify that we fill them in with zeros. By default it goes with NANs.
BD.restoreZeros(value=0)

# Export eigenvectors to files
print "Saving data to file"
f=open(f_out_path,'w')
for i in range (len(BD.dataDict['heatmap'])):
	f.write("\t".join(map(str,[genome_db_chrmLevel.idx2label[genome_db_chrmLevel.chrmIdxBinCont[i]],
															genome_db_chrmLevel.posBinCont[i],
															BD.EigDict['heatmap'][0][i].real]))+"\n")
f.close()

print "Finished!"