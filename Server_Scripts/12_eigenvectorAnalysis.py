"""
An example script to perform iterative correction of a raw heatmap and saving eigenvectors.
"""

import os
from mirnylib import genome
from mirnylib.h5dict import h5dict
from hiclib.binnedData import binnedData
import numpy as np

genome_name='mm9'

# When you do this, be sure that readChrms used to save heatmap matches
# readChrms that you define here!
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'Sp_full2'

domain_res=1000000

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5-raw'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'
f_out_path=base_folder+base_filename+"_"+str(domain_res/1000)+'.eig'

genome_fai_filepath='../fasta/'+genome_name+'/'+genome_name+'.fai'


NumEigenvectors = 1  # number of eigenvectors to compute

# Read resolution from one of the datasets
resolution = int(domain_res)

# Define the binnedData object, load data
BD = binnedData(resolution, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')

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
for i in range (len(BD.dataDict['HindIII_GM_1'][0])):
	strToWrite = ""
	curChrmIdx = genome_db.chrmIdxBinCont[i]
	if curChrmIdx == 0: 
		curRelativeBinNumb = i 
	else: 
		curRelativeBinNumb = i-genome_db.chrmLensBin[0:curChrmIdx].sum()
	
	if curChrmIdx < 19:
		chr_to_write='chr'+str(curChrmIdx+1)
	else:
		chr_to_write='chrX'
		
	strToWrite += chr_to_write+"\t"+str(genome_db.posBinCont[i])+"\t"+str(genome_db.posBinCont[i]+genome_db.binSizesBp[curChrmIdx][curRelativeBinNumb])
	strToWrite += "\t"+str(BD.EigDict['HindIII_GM_1'][0][i])+"\n"
	f.write (strToWrite)
f.close()

print "Finished!"
