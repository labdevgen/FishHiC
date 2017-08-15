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

from mirnylib.numutils import PCA, EIG, correct, \
    ultracorrectSymmetricWithVector, isInteger, \
    observedOverExpected, ultracorrect, adaptiveSmoothing, \
    removeDiagonals

genome_name='mm9'
genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'ESC_full'

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'
raw_heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'
figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb-Xist.png'

genome_fai_filepath='../../fasta/'+genome_name+'/'+genome_name+'.fai'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')
BD_raw=binnedData.binnedData(domain_res, genome_db)
BD_raw.simpleLoad(heatmap_filepath, 'HindIII_GM_1')


q=BD.dataDict['HindIII_GM_1']
q_raw=BD_raw.dataDict['HindIII_GM_1']

X_values=[]
Y_values=[]
Y_errors=[]

binnumber=100
dist=-1

start=sum(genome_db.chrmLensBin[0:19])
end=len(q)

print start
print end

for i in range(start,end):
	X_values.append(i-start+1)
	if abs(i-start-binnumber) <= dist:
		Y_values.append(-1.0)
		Y_errors.append(-1.0)
		continue
	else:
		Y_values.append(q[start+binnumber][i])
		Y_errors.append(1.0/(q_raw[start+binnumber][i]**0.5))

plt.subplots_adjust(bottom=0.15)
plt.plot(X_values,Y_values,"b-")

print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()

plt.clf()
plt.subplots_adjust(bottom=0.15)
plt.plot(X_values,Y_errors,"b-")

print "Saving figure "+figure_path+".errors.png"
f = open(figure_path+".errors.png", "wb")
plt.savefig(figure_path+".errors.png",dpi=600)
f.close()

