from mirnylib import genome
from hiclib import fragmentHiC
genome_name='mm9'
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
base_filename = 'SRR443884_mESC_2'
base_folder='/mnt/storage/home/vsfishman/HiC/data/'
fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'

#fragments = fragmentHiC.HiCdataset(
#     filename=fragment_dataset_filename,
#     genome=genome_db,
#     maximumMoleculeLength=500,
#     mode='r')

#fragments.CountDuplicateFragments()

print "Crating new fragments dataset\n"
fragments = fragmentHiC.HiCdataset(
     filename=fragment_dataset_filename,
     genome=genome_db,
     maximumMoleculeLength=500,
     mode='r')
# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
# at this stage, with maximumMoleculeLength specified at the initiation of the 
# object.
fragments.rebuildFragments()
print len(fragments.ufragments)