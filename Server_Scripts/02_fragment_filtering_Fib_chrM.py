from mirnylib import genome
from hiclib import fragmentHiC

#Define files and directories

base_filename = 'Fib_full2_chrM'
base_folder='/mnt/storage/home/vsfishman/HiC/data/'

fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'
heatmap_filepath=base_folder+'heatmap-res-'+str('fragment_')+'KB_'+base_filename+'.hdf5'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'

genome_name='mm9'
genome_folder='/mnt/storage/home/vsfishman/HiC/fasta/'
genome_db = genome.Genome(genome_folder+genome_name, readChrms=['M'])

metadata_file=base_folder+base_filename+"fragments_statistic.txt"

#genome_fai_filepath=genome_folder+genome_name+'/'+genome_name+'.fai'

# Create a HiCdataset object.
fragments = fragmentHiC.HiCdataset(
    filename=fragment_dataset_filename,
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='w')

# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
# at this stage, with maximumMoleculeLength specified at the initiation of the
# object.
fragments.parseInputData(
    dictLike=maped_reads_filepath)

fragments.setRfragAbsIdxs("HindIII")
#fragments.filterRsiteStart(offset=5)
#fragments.filterDuplicates()
#fragments.filterLarge()
#fragments.filterExtreme(cutH=0.005, cutL=0)

fragments.saveHeatmap(heatmap_filepath, resolution='fragment')

fragments.printMetadata(saveTo=metadata_file)