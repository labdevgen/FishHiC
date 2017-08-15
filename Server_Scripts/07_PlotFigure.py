import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lib_matrix_operations import *
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

#genome_name='mm9'
#genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
genome_name='mm9'
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=["#","X","Y"])

domain_res=200000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
#base_filename = 'SRR443884_mESC_2'
base_filename = 'Sp_full2_Y'
#base_filename = 'LA'


raw="-raw"
#raw=""

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'+raw
#heatmap_filepath=base_folder+'heatmap-res-'+str('fragment_')+'KB_'+base_filename+'.hdf5'+raw

#maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'

figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb'+raw+'.png'
#figure_path=base_folder+base_filename+"_"+str('fragment_')+'kb'+raw+'.eps'


BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')
q=np.log(BD.dataDict['HindIII_GM_1'])
#if domain_res < 1000000:
#	print "Matrix is too big, have to resize it"
#	q=resize_matrix(q,1000000/domain_res,np.max)
print "saving ",figure_path
plotting.plot_matrix(q)
plt.subplots_adjust(bottom=0.15)
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()
