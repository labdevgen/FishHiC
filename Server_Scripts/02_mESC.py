from mirnylib import genome
from hiclib import fragmentHiC 

genome_name='mm10'


# Create a HiCdataset object.
genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])
fragments = fragmentHiC.HiCdataset(
    filename='/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/fragment_dataset_mESC.hdf5',
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w')

# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
# at this stage, with maximumMoleculeLength specified at the initiation of the 
# object.
fragments.parseInputData(
    dictLike='/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/mapped_reads_mESC.hdf5')

fragments.filterRsiteStart(offset=5)
fragments.filterDuplicates()
fragments.filterLarge()
fragments.filterExtreme(cutH=0.005, cutL=0)

fragments.saveHeatmap('/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/heatmap-res-1M_mESC.hdf5', resolution=1000000)
