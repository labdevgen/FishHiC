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

domain_res=100000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
#base_filename1 = 'Fib_full2'
base_filename1 = 'Fib_full2'
#base_filename2 = 'Sp_full2'
base_filename2 = 'Sp_full2'
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

figure_path = "/mnt/storage/home/vsfishman/HiC/pics/"

chrms = range(genome_db1.chrmCount) #numner of chrms in genome
st = genome_db2.chrmStartsBinCont #array of numbers of chrms start (in bins)
end = genome_db2.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)

domains_file = "/mnt/storage/home/vsfishman/HiC/My_Scripts/nested_domains.txt"

big_domains = np.genfromtxt(domains_file,skip_header=1,filling_values="0.0",autostrip=True,usecols=(0,1,2),dtype=None)

domains_coverage = []

if not genome_db1.hasEnzyme():
	genome_db1.setEnzyme("HindIII")


for ind,chr in enumerate(big_domains['f0']):
	chr = genome_db1.label2idx[chr.split("chr")[-1]]

	ch_start = big_domains['f1'][ind]
	ch_end = big_domains['f2'][ind]
	 
	ch_start = ch_start / domain_res
	ch_end = ch_end / domain_res
	i=chr
	q1_nonzero_c=np.count_nonzero(q1[st[i]:end[i],st[i]:end[i]])
	q2_nonzero_c=np.count_nonzero(q2[st[i]:end[i],st[i]:end[i]])
	
	domains_coverage.append([float(q2_nonzero_c)/float(q1_nonzero_c),float(q1_nonzero_c)/(ch_end-ch_start),float(ch_end-ch_start),chr,ch_start,ch_end])

domains_coverage=np.array(domains_coverage)

print domains_coverage
domains_coverage=domains_coverage.view(np.recarray).view(dtype=[('ratio', '<f8'),('length_ratio', '<f8'), ('length', '<f8'), ('chr', '<f8'), ('st', '<f8'), ('end', '<f8')])
domains_coverage_length = np.array(domains_coverage)
domains_coverage.sort(order=['ratio'], axis=0)
domains_coverage_length.sort(order=['length'], axis=0)
print "After sort:"
print domains_coverage

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
	plotting.plot_matrix(np.log(q2_0),cmap='OrRd')
	plt.subplots_adjust(bottom=0.15)
	fp1=figure_path+str(domain_res)+"_1000KB_chr"+str(i+1)+".png"
	print "Saving figure "+fp1
	f = open(fp1, "wb")
	plt.savefig(fp1,dpi=900)
	f.close()
	plt.clf()

def plot_one_chr_fragment(matrix,figure_path,base_filename,chr_numb,ch_start,ch_end,domain_st,domain_end): #chr_numb is zero-based
	print "Plotting picture"
	ch_start = ch_start / domain_res
	ch_end = ch_end / domain_res
	i=chr_numb
	q2_0=matrix[st[i]:end[i],st[i]:end[i]]
	q2_0=q2_0[ch_start:ch_end,ch_start:ch_end]
	plotting.plot_rotated_matrix(np.log(q2_0),45,cmap='OrRd')
	print ch_start,ch_end,domain_st,domain_end
	print domain_st-ch_start,domain_end-ch_start
	plt.axvline(x=domain_st-ch_start)
	plt.axvline(x=domain_end-ch_start)
	plt.subplots_adjust(bottom=0.15)
	fp1=figure_path+str(domain_res)+"_chr"+str(i+1)+"_bin_"+str(ch_start)+"_to_"+str(ch_end)+"_"+base_filename+".rotated.png"
	print "Saving figure "+fp1
	f = open(fp1, "wb")
	plt.savefig(fp1,dpi=900)
	f.close()
	plt.clf()

#plot_one_chr_pict(q2,figure_path+base_filename2,18)

domains_coverage=domains_coverage_length
borders = 2
for i in xrange(len(domains_coverage['chr'])-20,len(domains_coverage['chr'])): #plot figures for top20 domains
	s = max (0,domains_coverage['st'][i][0]-borders)
	e = min (genome_db2.chrmLensBin[domains_coverage['chr'][i][0]]-1,domains_coverage['end'][i][0]+borders)
	plot_one_chr_fragment(q1,figure_path,base_filename1,domains_coverage['chr'][i][0],s*domain_res,e*domain_res,domains_coverage['st'][i][0],domains_coverage['end'][i][0])
	plot_one_chr_fragment(q2,figure_path,base_filename2,domains_coverage['chr'][i][0],s*domain_res,e*domain_res,s+domains_coverage['st'][i][0],domains_coverage['end'][i][0])