print "Starting programm"
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lib_matrix_operations import *
import scipy.stats
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC
import math

genome_name='mm9'

genome_db1 = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
genome_db2 = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
#base_filename1 = 'Fib_full2'
base_filename1 = 'Sp_full2'
#base_filename2 = 'Sp_full2'
base_filename2 = 'Fib_full2'
#base_filename2 = 'Fib_full_compressed_as_ESC_full'

#IMPORTANT: use iter-corrected heatmaps here. Otherwise, take care about adjustment of total reads number when calculating mask_hugeDifference
heatmap_filepath1=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename1+'.hdf5'
heatmap_filepath2=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'.hdf5'

print "Loading file "+heatmap_filepath1
BD1 = binnedData.binnedData(domain_res, genome_db1)
BD1.simpleLoad(heatmap_filepath1, 'HindIII_GM_1')

print "Loading file "+heatmap_filepath2
BD2 = binnedData.binnedData(domain_res, genome_db2)
BD2.simpleLoad(heatmap_filepath2, 'HindIII_GM_1')

q1=BD1.dataDict['HindIII_GM_1']
q2=BD2.dataDict['HindIII_GM_1']

def calcSpearm(q1,q2,BD1):
	mask = (np.sum(q1, axis=0) > 0) * (np.sum(q2, axis=0) > 0)
	validMask = mask[:, None] * mask[None, :]
	transmask = BD1.chromosomeIndex[:, None] == BD1.chromosomeIndex[None, :]
	cormask = transmask * validMask
	c = scipy.stats.spearmanr(q1[cormask],q2[cormask])
	print c

#calcSpearm(q1,q2,BD1)

figure_path = "/mnt/storage/home/vsfishman/HiC/pics/"

chrms = range(genome_db1.chrmCount) #numner of chrms in genome
st = genome_db2.chrmStartsBinCont #array of numbers of chrms start (in bins)
end = genome_db2.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)

#print "Plotting error matrix"
#plotting.plot_matrix(np.log(q1))
#plt.subplots_adjust(bottom=0.15)
#fp1=figure_path+base_filename1+"_1MB.png"
#print "Saving figure "+fp1
#f = open(fp1, "wb")
#plt.savefig(fp1,dpi=600)
#f.close()
#plt.clf()

#print "Plotting error matrix"
#plotting.plot_matrix(np.log(q2))
#plt.subplots_adjust(bottom=0.15)
#fp1=figure_path+base_filename2+"_1MB.png"
#print "Saving figure "+fp1
#f = open(fp1, "wb")
#plt.savefig(fp1,dpi=600)
#f.close()
#plt.clf()

#i=18
#q1_0=q1[st[i]:end[i],st[i]:end[i]]
#plotting.plot_matrix(np.log(q1_0))
#plt.subplots_adjust(bottom=0.15)
#fp1=figure_path+base_filename1+"_1MB_chr"+str(i+1)+".png"
#print "Saving figure "+fp1
#f = open(fp1, "wb")
#plt.savefig(fp1,dpi=600)
#f.close()
#plt.clf()

def plot_one_chr_pict(matrix,figure_path,chr_numb): #chr_numb is zero-based
	print "Plotting picture"
	i=chr_numb
	q2_0=matrix[st[i]:end[i],st[i]:end[i]]
	plotting.plot_matrix(np.log(q2_0))#,cmap='OrRd')
	plt.subplots_adjust(bottom=0.15)
	fp1=figure_path+str(domain_res)+"_1000KB_chr"+str(i+1)+".png"
	print "Saving figure "+fp1
	f = open(fp1, "wb")
	plt.savefig(fp1,dpi=900)
	f.close()
	plt.clf()

def plot_one_chr_fragment(matrix,figure_path,chr_numb,ch_start,ch_end): #chr_numb is zero-based
	print "Plotting picture"
	ch_start = ch_start / domain_res
	ch_end = ch_end / domain_res
	i=chr_numb
	q2_0=matrix[st[i]:end[i],st[i]:end[i]]
	q2_0=q2_0[ch_start:ch_end,ch_start:ch_end]
	plotting.plot_matrix(np.log(q2_0),cmap='OrRd')
	plt.subplots_adjust(bottom=0.15)
	fp1=figure_path+str(domain_res)+"_40KB_chr"+str(i+1)+"_bin_"+str(ch_start)+"_to_"+str(ch_end)+".png"
	print "Saving figure "+fp1
	f = open(fp1, "wb")
	plt.savefig(fp1,dpi=900)
	f.close()
	plt.clf()
for ch in [1,2]:
	plot_one_chr_pict(q1,figure_path+base_filename1,ch)
	plot_one_chr_pict(q2,figure_path+base_filename2,ch)