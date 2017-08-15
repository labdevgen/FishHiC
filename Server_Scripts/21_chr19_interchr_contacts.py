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
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'Fib_full2'

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'
chr_numb = 18
figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+"kb-intrarchr-chr"+str(chr_numb+1)+".png"

genome_fai_filepath='../fasta/'+genome_name+'/'+genome_name+'.fai'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')
q=BD.dataDict['HindIII_GM_1']

print "Calculating"

chr_start_bin=sum(genome_db.chrmLensBin[0:chr_numb])
chr_end_bin=sum(genome_db.chrmLensBin[0:chr_numb+1])
step = 5
l = (chr_end_bin - chr_start_bin)
res = np.zeros(shape=(l/step  + 1,genome_db.chrmCount))

for i in xrange(0,l,step):
	for j in range(genome_db.chrmCount):
		if j == chr_numb: continue	#remove contacts with cromosome itself

		chr_j_start_bin=sum(genome_db.chrmLensBin[0:j])
		chr_j_end_bin=sum(genome_db.chrmLensBin[0:j+1])

		res[(i/step),j] = np.sum(q[(i+chr_start_bin):(i+chr_start_bin+step),chr_j_start_bin:chr_j_end_bin])

total=np.sum(res)

res_sums = [sum(res[:,j]) for j in range(genome_db.chrmCount)]

res_probablities = np.zeros_like(res)

for i in xrange(0,l,step):
	i=(i/step)
	for j in xrange(genome_db.chrmCount):
		if j == chr_numb: 
			res_probablities[i,j] = None
			continue	#remove contacts with cromosome itself

		s1=float(np.sum(res[i,:]))
		s2=float(np.sum(res[:,j]))
		p1=(s1/total)*(s2/(total-s1))
		p2=(s2/total)*(s1/(total-s2))

		znam = (p1+p2)*(total/2.0)
		if znam == 0:
			res_probablities[i,j] = None
		else:
			res_probablities[i,j]  = res[i,j]/znam

#res_probablities=np.log2(res_probablities)

mi=np.min(res_probablities) #0.8 #Min and Max values for pictures
ma=np.max(res_probablities) #1.4 #To make colors same
if (np.max(res_probablities) > ma) or (np.min(res_probablities) < mi):
	print "Current max and min are ",np.max(res_probablities),np.min(res_probablities)
	raise Exception("Array values out of "+str(mi)+"\t"+str(ma))

print "Plotting contact matrix"
#plotting.plot_matrix(res_probablities,vmin=mi, vmax=ma)
plotting.plot_matrix(res_probablities)
plt.subplots_adjust(bottom=0.15)
print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()
