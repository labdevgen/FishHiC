import datetime

now = datetime.datetime.now()

print "Starting domains_search_Dixon "
print str(now)
print "\n"

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

genome_name='mm10'
genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res = 100000

#Define files and directories
base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'SRR443883_mESC_1'
HMM_file_path='/mnt/storage/home/vsfishman/HiC/DI/domaincall_software/'

fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'
heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'
figure_path=base_folder+base_filename+str(domain_res/1000)+'kb.png'

genome_fai_filepath='../../fasta/'+genome_name+'/genome_1.fai'

def executeBashCommand(bashCommand):
  print "executing command %s \n" % (bashCommand)
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  print p.communicate()[0]

def doDixonSearch():
  bashCommand ="./DI_from_matrix_minja.pl "+heatmap_filepath+".matrix "+str(domain_res)+" "+str(domain_res*50)+" "+genome_fai_filepath+" "+heatmap_filepath+".matrix.di"
  executeBashCommand(bashCommand)
  print "\n Done!\n"

def doDixon_PostProc():
  print "----Starting dixon postproc---\n"
  print "Generating new HMM_calls file ("+HMM_file_path+'HMM_calls_'+base_filename+".m)\n"
  m_inp_file = open('HMM_calls.m','r')
  m_out_file = open('HMM_calls_'+base_filename+'.m','w')
  m_content = m_inp_file.read()
  m_content = m_content.replace('i_filename',heatmap_filepath+".matrix.di")
  m_content = m_content.replace("addpath(genpath('","addpath(genpath('"+HMM_file_path)

  m_out_file.write(m_content)
  m_inp_file.close()
  m_out_file.close()
#  executeBashCommand("bash ~/modules.sh'")
#  executeBashCommand("module load matlab")
  bashCommand ="matlab nodisplay -nodesktop -nojvm -nosplash -r HMM_calls_"+base_filename
  executeBashCommand(bashCommand)
  bashCommand="perl "+HMM_file_path+"/perl_scripts/file_ends_cleaner.pl "+heatmap_filepath+".matrix.di.m"+heatmap_filepath+".matrix.di | perl "+"/perl_scripts/converter_7col.pl > "+heatmap_filepath+".matrix.di.m.hmm_7colfile"  
  executeBashCommand(bashCommand)
  print "\n Done!\n"
  


###################Step 1 - create heatmap with defined resolution################

if (os.path.isfile(heatmap_filepath+".matrix.di")):
   print ".di file found, assuming it is already correct\n"
   doDixon_PostProc()
   exit()

if (os.path.isfile(heatmap_filepath+".matrix")):
   print "Matrix file found, assuming it is already correct\n"
   doDixonSearch()
   exit()

if (os.path.isfile(heatmap_filepath)):
   print "heatmap %s aready exists, skipping...\n" % (heatmap_filepath)
elif (os.path.isfile(fragment_dataset_filename)):
   print "fragments dataset found, assuming it is already correct\n"	
   # Create a HiCdataset object.
   fragments = fragmentHiC.HiCdataset(
     					 filename=fragment_dataset_filename,
    					  genome=genome_db,
    					  maximumMoleculeLength=500,
    					  mode='r')	
   fragments.saveHeatmap(heatmap_filepath, domain_res)
else:
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

   fragments.filterRsiteStart(offset=5)
   fragments.filterDuplicates()
   fragments.filterLarge()
   fragments.filterExtreme(cutH=0.005, cutL=0)

   fragments.saveHeatmap(heatmap_filepath, domain_res)

###################Step 2 - create Matrix################

# Read resolution from the dataset.
raw_heatmap = h5dict.h5dict(heatmap_filepath, mode='r') 
resolution = int(raw_heatmap['resolution'])

####### Set resolution for genome
#genome_db.setResolution(resolution)

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
#BD.iterativeCorrectWithoutSS()

####### Perform fake cis.
#BD.fakeCis()

# Plot the heatmap directly.
#plotting.plot_matrix(np.log(BD.dataDict['HindIII_GM_1']))

#plt.show()
#plt.savefig(figure_path)

print "Writing file %s with the information \n"% (heatmap_filepath+'.matrix')
print "ChrIndex \t StartBin(Nucleotyde) \t EndBin(Nucl) \t Values\n"

f=open(heatmap_filepath+'.matrix','w')

for i in range (len(BD.dataDict['HindIII_GM_1'])):
    strToWrite = ""
    curChrmIdx = genome_db.chrmIdxBinCont[i]
    if curChrmIdx == 0: 
	curRelativeBinNumb = i 
    else: 
	curRelativeBinNumb = i-genome_db.chrmLensBin[0:curChrmIdx].sum()

    strToWrite += str(curChrmIdx)+"\t"+str(genome_db.posBinCont[i])+"\t"+str(genome_db.posBinCont[i]+genome_db.binSizesBp[curChrmIdx][curRelativeBinNumb])
    for j in range (len(BD.dataDict['HindIII_GM_1'])): 
       strToWrite += "\t"+str(BD.dataDict['HindIII_GM_1'][i][j])
    strToWrite += "\n"
    f.write (strToWrite)

f.close()

doDixonSearch()
