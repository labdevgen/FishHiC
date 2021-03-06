from mirnylib import genome
from hiclib import fragmentHiC 

# Create a HiCdataset object.
genome_db = genome.Genome('../../fasta/mm10', readChrms=['#', 'X'])
fragments = fragmentHiC.HiCdataset(
    filename='../../data/serov/fragment_dataset_Fib.hdf5',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w')

# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
# at this stage, with maximumMoleculeLength specified at the initiation of the 
# object.
fragments.parseInputData(
    dictLike='../../data/serov/mapped_reads_Fib.hdf5')

fragments.filterRsiteStart(offset=5)
fragments.filterDuplicates()
fragments.filterLarge()
fragments.filterExtreme(cutH=0.005, cutL=0)

fragments.saveHeatmap('../../data/serov/heatmap-res-1M_Fib.hdf5', resolution=1000000)
