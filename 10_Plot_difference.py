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

def normalize(q):
	s=0
	i=0
	while (s==0):
		i +=1
		s=sum(q[i])
	return q/s
	
genome_name='mm9'
genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename1 = 'Fib_full'
base_filename2 = 'Sp_full'

heatmap_filepath1=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename1+'.hdf5'
heatmap_filepath2=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'.hdf5'

figure_path=base_folder+"difference"+base_filename1+'-'+base_filename2+str(domain_res/1000)+'kb.png'

print "Creating objects"
BD1 = binnedData.binnedData(domain_res, genome_db)
BD2 = binnedData.binnedData(domain_res, genome_db)

print "Loading data"
BD1.simpleLoad(heatmap_filepath1, 'HindIII_GM_1')
BD2.simpleLoad(heatmap_filepath2, 'HindIII_GM_1')


print "Normalizing data"
q1=normalize(BD1.dataDict['HindIII_GM_1'])
q2=normalize(BD2.dataDict['HindIII_GM_1'])

print "Plotting data"
plotting.plot_matrix(np.log(q1/q2))

print "Saving figure"
plt.subplots_adjust(bottom=0.15)
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()
print "Done!"