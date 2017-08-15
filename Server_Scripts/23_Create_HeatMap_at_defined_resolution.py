import datetime

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--res")
parser.add_argument("--base_filename")
args = parser.parse_args()

domain_res=args.res
base_filename=str(args.base_filename).strip()

if (domain_res==None) or (domain_res==""):
  domain_res=100000
if (base_filename==None):
  base_filename = ''

domain_res=int(domain_res)

print "Domain resolution = "+str(domain_res)


#Define files and directories
base_folder='/mnt/storage/home/vsfishman/HiC/data/'

fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'
heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'

genome_name='mm9'
genome_folder='/mnt/storage/home/vsfishman/HiC/fasta/'
genome_db = genome.Genome(genome_folder+genome_name, readChrms=['#', 'X'])
genome_fai_filepath=genome_folder+genome_name+'/'+genome_name+'.fai'

def CorrectHeatMap():
  # Read resolution from the dataset.
  print "Loading raw heatmap\n"
  raw_heatmap = h5dict.h5dict(heatmap_filepath+'-raw', mode='r') 
  resolution = int(raw_heatmap['resolution'])

  ####### Set resolution for genome
  #genome_db.setResolution(resolution)

  # Create a binnedData object, load the data.
  BD = binnedData.binnedData(resolution, genome_db)
  BD.simpleLoad(heatmap_filepath+'-raw', 'HindIII_GM_1')
  # Remove the contacts between loci located within the same bin.
  BD.removeDiagonal()
  # Remove bins with less than half of a bin sequenced.
  BD.removeBySequencedCount(0.5)
  # Remove 1% of regions with low coverage.
  BD.removePoorRegions(cutoff=1)
  # Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
  BD.truncTrans(high=0.0005)
  # Perform iterative correction.
  BD.iterativeCorrectWithoutSS(force=True)
  # Save the iteratively corrected heatmap.
  BD.export('HindIII_GM_1', heatmap_filepath)

  ####### Perform fake cis. 
  #BD.fakeCis()

  # Plot the heatmap directly.
  #plotting.plot_matrix(np.log(BD.dataDict['HindIII_GM_1']))

  #plt.show()
  #plt.savefig(figure_path)

fragments = fragmentHiC.HiCdataset(
     filename=fragment_dataset_filename,
     genome=genome_db,
     maximumMoleculeLength=500,
     mode='r')

fragments.saveHeatmap(heatmap_filepath+'-raw', domain_res)
CorrectHeatMap()