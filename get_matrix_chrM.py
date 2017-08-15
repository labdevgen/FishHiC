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
from mirnylib.h5dict import h5dict


genome_name='mm9'
#genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
#genome_name='hg19'
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['M'])

domain_res='fragment'

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
#base_filename = 'SRR443884_mESC_2'
base_filename = 'Fib_full2_chrM'

#raw="-raw"
raw=""

#heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'+raw
heatmap_filepath=base_folder+'heatmap-res-fragment_KB_'+base_filename+'.hdf5'+raw

#figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb'+raw+'.eps'
figure_path=base_folder+base_filename+'_fragment_kb'+raw+'.eps'


path = os.path.abspath(os.path.expanduser(heatmap_filepath))
alldata = h5dict(path , mode="r")

alldata["heatmap"]