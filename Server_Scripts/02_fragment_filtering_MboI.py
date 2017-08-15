from mirnylib import genome
from hiclib import fragmentHiC
import numpy as np

#Define files and directories

base_filename = 'Kalhor_Mitoch'
base_folder='/mnt/storage/home/vsfishman/HiC/ira/'

fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'
heatmap_filepath=base_folder+'heatmap-res-'+str('fragment_')+'KB_'+base_filename+'.hdf5'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'

genome_name='hg19'
genome_folder='/mnt/storage/home/vsfishman/HiC/fasta/'
genome_db = genome.Genome(genome_folder+genome_name, readChrms=['#','X','M'])

metadata_file=base_folder+base_filename+"fragments_statistic.txt"

#genome_fai_filepath=genome_folder+genome_name+'/'+genome_name+'.fai'

# Create a HiCdataset object.
fragments = fragmentHiC.HiCdataset(
    filename=fragment_dataset_filename,
    genome=genome_db,
    maximumMoleculeLength=500,
    mode='r')

# Load the parsed reads into the HiCdataset. The dangling-end filter is applied
# at this stage, with maximumMoleculeLength specified at the initiation of the
# object.
#fragments.parseInputData(
#    dictLike=maped_reads_filepath,enzymeToFillRsites="MboI")

#fragments.setRfragAbsIdxs("MboI")
#fragments.filterRsiteStart(offset=5)
#fragments.filterDuplicates()
#fragments.filterLarge()
#fragments.filterExtreme(cutH=0.005, cutL=0)

#fragments.saveHeatmap(heatmap_filepath, resolution='fragment')

#fragments.printMetadata(saveTo=metadata_file)

print "Generating rfragAbsIdxs"
fragments.genome.setEnzyme("MboI")

chrM_idx=fragments.genome.label2idx['M']
print chrM_idx
idxs = np.logical_and((fragments.chrms1 == chrM_idx),(fragments.chrms2 == chrM_idx))
print np.sum(idxs)

fragids1 = np.array(fragments.fragids1[idxs])
fragids2 = np.array(fragments.fragids2[idxs])

rfragAbsIdxs1 = fragments.genome.getRfragAbsIdxs(fragids1) #all rfgrags belonging to chrM
rfragAbsIdxs2 = fragments.genome.getRfragAbsIdxs(fragids2)

min_rfragAbsIdx = min(np.min(rfragAbsIdxs1),np.min(rfragAbsIdxs2))

print np.min(rfragAbsIdxs1)
print np.max(rfragAbsIdxs1)

print np.min(rfragAbsIdxs2)
print np.max(rfragAbsIdxs2)

print min_rfragAbsIdx

rfragAbsIdxs1 = rfragAbsIdxs1-min_rfragAbsIdx
rfragAbsIdxs2 = rfragAbsIdxs2-min_rfragAbsIdx

print np.min(rfragAbsIdxs1)
print np.max(rfragAbsIdxs1)

print np.min(rfragAbsIdxs2)
print np.max(rfragAbsIdxs2)


all_r_sites_on_M = fragments.genome.rsites[chrM_idx]

print "Calculating label array"
numBins = len(all_r_sites_on_M)
label = rfragAbsIdxs1
print label
label *= numBins
print label
label += rfragAbsIdxs2
print label

print "Counting label array"
minl = numBins**2
print minl
print max(label)
print label	
counts = np.bincount(label, minl)
print numBins
print len(label)
print len(counts)
print counts
print "Reshaping label array"
counts.shape = (numBins, numBins)
print "Calculating heatmap"
for i in xrange(len(counts)):
    counts[i, i:] += counts[i:, i]
    counts[i:, i] = counts[i, i:]

print "Exporting heatmap"
np.savetxt(base_folder+base_filename+'.chrM_textheatmap')
