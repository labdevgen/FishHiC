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

genome_name='mm9'
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
#base_filename = 'SRR443884_mESC_2'
base_filename = 'Sp_full2'

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5-raw'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'
figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb-errors.png'

genome_fai_filepath='../../fasta/'+genome_name+'/'+genome_name+'.fai'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')
q=BD.dataDict['HindIII_GM_1']
print "Elocating memory for error matrix"
q1=np.array(q)
print "Calculating error matrix"
s=0.0 #total summ of errors
n = 0 #number of mappable dots

chrms = range(genome_db.chrmCount) #numner of chrms in genome

def calc_one_chr(matrix,chr_numb,max_range=0): #chr_numb is zero-based
	print "Plotting picture"
	i=chr_numb
	st = genome_db.chrmStartsBinCont #array of numbers of chrms start (in bins)
	end = genome_db.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)
	q2=matrix[st[i]:end[i],st[i]:end[i]]
	s = 0
	n = 0
	for i in xrange(len(q2)):
		if max_range != 0:
			end = min (len(q2),i+max_range)
		for j in xrange(i+1,end):
			if q[i,j]==0: continue
			s += np.sqrt(q[i,j])/q[i,j]
			n += 1
	return (s,n)

for i in chrms:
	print "chr",i+1
	r = calc_one_chr(q,i,max_range=50)
	s += r[0]
	n += r[1]

print "Average error for interchromosomal contacts only = ", s/n

s=0.0 #total summ of errors
n = 0 #number of mappable dots

for i in range(len(q)):
	if i % (len(q) / 10) == 0:
		print (i / (len(q) / 10))*10,"%"
	for j in range(len(q[i])):
		if q[i,j]==0: 
			q1[i,j]=1
		else:
			q1[i,j]=1.0/np.sqrt(q[i,j])
			s += np.sqrt(q[i,j])/q[i,j]
			n += 1

print "Average error for all contacts only = ", s/n

if domain_res < 1000000:
	print "Matrix is too big, have to resize it"
	q1=resize_matrix(q1,1000000/domain_res,np.min)

print "Plotting error matrix"
plotting.plot_matrix(q1)
plt.subplots_adjust(bottom=0.15)
print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()