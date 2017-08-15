import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

genome_name='mm9'
genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'Sp_full'

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5-raw'

genome_fai_filepath='../fasta/'+genome_name+'/'+genome_name+'.fai'

BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')

sums=[]

for j in range(len(BD.dataDict['HindIII_GM_1'])):
  sums.append(float(sum(BD.dataDict['HindIII_GM_1'][i])))
maxs = max(sums)
sums = np.array(sums)/maxs

for j in range(len(BD.dataDict['HindIII_GM_1'])):
  
  sums.append(float(sum(BD.dataDict['HindIII_GM_1'][i])))


plotting.plot_matrix(q)
plt.savefig(figure_path)