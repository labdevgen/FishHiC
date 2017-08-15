import matplotlib.pyplot as plt
import numpy as np

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData

genome_name='mm10'
genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])
heatmap_filepath='/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/heatmap-res-1M_mESC.hdf5'


# Read resolution from the dataset.
raw_heatmap = h5dict.h5dict(heatmap_filepath, mode='r') 
resolution = int(raw_heatmap['resolution'])

# Create a binnedData object, load the data.
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')

# Remove the contacts between loci located within the same bin.
BD.removeDiagonal()

# Remove bins with less than half of a bin sequenced.
BD.removeBySequencedCount(0.5)

# Remove 1% of regions with low coverage.
BD.removePoorRegions(cutoff=1)

# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
BD.truncTrans(high=0.0005)

# Perform iterative correction.
BD.iterativeCorrectWithoutSS()

# Save the iteratively corrected heatmap.
BD.export('HindIII_GM_1', heatmap_filepath+'_corrected')

# Plot the heatmap directly.
plotting.plot_matrix(np.log(BD.dataDict['HindIII_GM_1']))

#plt.show()
plt.savefig('mESC.png')
