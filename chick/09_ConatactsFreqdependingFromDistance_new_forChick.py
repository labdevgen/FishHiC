import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import math
import numpy as np
import scipy

from mirnylib import genome
from mirnylib import h5dict
from hiclib import binnedData


colors=['b','g','r','c','m','y','k','w']
markers=['o','v','s','p','*','+','x','D']
chrms_in_bins=[] #Item number i of this array is chromosome number for bin number i
genome_db=None
resolution=0

filenames={}
filenames['/mnt/storage/home/vsfishman//HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII-100k.hm']='Blood100Kb'
filenames['/mnt/storage/home/vsfishman//HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII-100k.hm']='ChEF100Kb'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_ESC_full.hdf5']='ES_full'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_SRR443884_mESC_2.hdf5']='ES_2'
filenames_sums={}


use_weights=True
split_by_chromosomes=True

def format_plot_axes():
	for i in [1]:
		plt.figure(i)
		ax=plt.gca()
		plt.setp(ax.get_xticklabels(), rotation='vertical')
		new_ticks=ax.get_xticks()
		additional_ticks=np.arange(new_ticks[0]+new_ticks[1]/5,new_ticks[1],new_ticks[1]/5)
		new_ticks=sorted(additional_ticks.tolist()+new_ticks.tolist())
		plt.xticks(new_ticks,[str(x/1000)+"KB" if x < 1000000 else str(x/1000000)+"MB" for x in new_ticks], rotation='vertical')
		ax.set_xlim(0,max(new_ticks)+1)
		ax.yaxis.grid(b=True, which='both', color='gray', linestyle='dotted',linewidth=0.5)
		ax.xaxis.grid(b=True, which='both', color='gray', linestyle='dotted',linewidth=0.5)
		if use_weights:
			plt.ylabel('% of contacts')
		else:
			plt.ylabel('Number of contacts')
		plt.xlabel('Distance from current bin')
		plt.title(str(filenames.values())+" weights="+str(use_weights))	
		plt.subplots_adjust(bottom=0.15)
#		ax.set_yscale('log',base=2.718281828)
#		ax.set_xscale('log',base=2.718281828)

diags={}

for i in filenames.keys():  
	print "Reading file "+i
#	if (i.split('.')[-1]=='hdf5'):
	if True:
		if (resolution==0): #if we do not know resolution
			raw_heatmap = h5dict.h5dict(i, mode='r') #open heatmap
			resolution = int(raw_heatmap['resolution']) #get the resolution
			del raw_heatmap #close heatmap
		if (genome_db==None): #if we have not initilaize genome before
			genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")
		BD = binnedData.binnedData(resolution, genome_db) #now we can initialyze heatmap with defined resolution and genome
		BD.simpleLoad(i, 'heatmap')
		number_of_bins=len(BD.dataDict['heatmap'])
		
		diags[i]=np.zeros(max(genome_db.chrmLensBin))
		for chr in xrange(genome_db.chrmCount):
			for j in xrange(genome_db.chrmLensBin[chr]):
				cur_chr_matrix = BD.dataDict['heatmap'][genome_db.chrmStartsBinCont[chr]:genome_db.chrmEndsBinCont[chr],genome_db.chrmStartsBinCont[chr]:genome_db.chrmEndsBinCont[chr]]
				diags[i][j] += sum(np.diag(cur_chr_matrix,j))*2
		print np.sum(diags[i])
		diags[i] = (diags[i]/np.sum(diags[i]))*100.0
		print (diags[i][0:10])/100.0
		del BD		


def contact_freq_total():
	X_min = -1
	for i,val in enumerate(diags[diags.keys()[0]]):
		if ((val != 0) and (X_min == -1)):
			X_min = i
		if ((X_min != -1) and (val <= diags[diags.keys()[0]][X_min]*0.001)):
			X_max = i
			break
			
	
	color_numb = 0
	print X_max, X_min
	for i in sorted(filenames.keys()):
			plt.figure(1)
			line_format=colors[color_numb]+markers[color_numb]+'-'
			plt.plot([x*resolution for x in range(X_min,X_max)],diags[i][X_min:X_max],
					line_format,markersize=1,
					color=colors[color_numb],label=filenames[i])
			color_numb +=1
			if (color_numb >= len(colors)-1 or color_numb>=len(markers)-1):
				color_numb = 0 
	
	format_plot_axes()
	plt.figure(1)
	plt.legend(loc="upper right",bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	filename="Distribution_result_total_"+str(resolution/1000)+"KB"
	if use_weights:
		filename += "_Weighted"
	filename+=".png"
	f = open(filename, "wb")
	plt.savefig(filename,format='png',dpi=300)
	f.close()
	plt.clf()

def contact_freq_by_chromosomes():
	X_values_list=range(1,number_of_bins+1)
	chromosomes_list=[0,18,19]
	for i in sorted(filenames.keys()):
		color_numb=0
#		for j in contact_freq[i].keys():
		for j in chromosomes_list:			
			plt.figure(1)
			contact_length_90=len(contact_freq[i][j]["In"])
			line_format=colors[color_numb]+markers[color_numb]+'-'
			plt.plot(X_values_list[1:contact_length_90-1],contact_freq[i][j]["In"][1:contact_length_90-1],
					line_format,markersize=1,
					color=colors[color_numb],label=filenames[i]+"_chr"+str(j+1))
			color_numb +=1
			if (color_numb >= len(colors)-1 or color_numb>=len(markers)-1):
				color_numb = 0 
		format_plot_axes()
		plt.figure(1)
		plt.legend(loc="upper right",bbox_to_anchor=(1.05, 1), borderaxespad=0.)
		filename="Distribution_result_byCHR_"+str(resolution/1000)+"KB_"+filenames[i]
		if use_weights:
			filename += "_Weighted"
		filename+=".png"
		f = open(filename, "wb")
		plt.savefig(filename,format='png',dpi=300)
		f.close()
		plt.clf()

def contact_sums_total():
	plt.clf()
	plt.figure(1)
	X_min = -1
	for i,val in enumerate(diags[diags.keys()[0]]):
		if ((val != 0) and (X_min == -1)):
			X_min = i
		if ((X_min != -1) and (val <= diags[diags.keys()[0]][X_min]*0.001)):
			X_max = i
			break
	
	color_numb = 0
	X_min = 0
	X_max = 200
	plt.gca().set_xlim(0,X_max*resolution+1)
	plt.gca().set_xticks(range(0,X_max*resolution,X_max*resolution/10))
	
	for i in sorted(filenames.keys()):
			plt.figure(1)
			line_format=colors[color_numb]+markers[color_numb]+'-'
			plt.plot([x*resolution for x in range(X_min,X_max)],np.cumsum(diags[i][X_min:X_max]),
					line_format,markersize=1,
					color=colors[color_numb],label=filenames[i])
			print plt.gca().get_xticks()
			color_numb +=1
			if (color_numb >= len(colors)-1 or color_numb>=len(markers)-1):
				color_numb = 0 

	format_plot_axes()
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
	filename="Distribution_result_sums_total_"+str(resolution/1000)+"KB"
	if use_weights:
		filename += "_Weighted"
	filename+=".png"
	f = open(filename, "wb")
	plt.savefig(filename,format='png',dpi=300)
	f.close()
	plt.clf()


def contact_sums_by_chromosomes():
	X_values_list=range(1,number_of_bins+1)
	chromosomes_list=[0,18,19]
	contact_length_90=min(genome_db.chrmLensBin)+1
	
	for i in sorted(filenames_sums.keys()):
		color_numb=0
#		for j in contact_freq[i].keys():
		for j in chromosomes_list:			
			plt.figure(1)
			contact_length_90=genome_db.chrmLensBin[j]+1
			line_format=colors[color_numb]+markers[color_numb]+'-'
			plt.plot(X_values_list[1:contact_length_90-1],contact_freq_sums[i][j]["In"][1:contact_length_90-1],
					line_format,markersize=1,
					color=colors[color_numb],label=filenames_sums[i]+"_chr"+str(j+1))
			color_numb +=1
			if (color_numb >= len(colors)-1 or color_numb>=len(markers)-1):
				color_numb = 0 
		format_plot_axes()
		plt.figure(1)
		plt.legend(loc="center right",bbox_to_anchor=(1.05, 1), borderaxespad=0.)
		filename="Distribution_result_byCHR_"+str(resolution/1000)+"KB_"+filenames_sums[i]
		if use_weights:
			filename += "_Weighted"
		filename+=".png"
		f = open(filename, "wb")
		plt.savefig(filename,format='png',dpi=300)
		f.close()
		plt.clf()


#Printing contact freq figure
#contact_freq_by_chromosomes()
contact_freq_total()
#contact_sums_by_chromosomes()
contact_sums_total()

#color_numb=0

#for i in sorted(filenames_sums.keys()):
#	plt.figure(2)
#	line_format=colors[color_numb]+'-'
#	plt.plot(X_values_list[1:],contact_freq_sums[i][1:],
#			line_format,markersize=1,
#			color=colors[color_numb],label=filenames_sums[i])
#	color_numb +=1
#	if (color_numb >= len(colors)-1 or color_numb>=len(markers)-1):
#		color_numb = 0 


#plt.figure(2)
#plt.legend(loc=1)
#f = open('Distribution_result_sum-100KB-noWeight.png', "wb")
#plt.savefig('Distribution_result_sum-100KB.png',format='png',dpi=300)	
#f.close()
