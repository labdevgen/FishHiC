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
genome_name='mm9'
genome_db=None
resolution=0

filenames={}
filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_Fib_full.hdf5']='FibMatr1MB'
filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_Sp_full.hdf5']='SpMatr1MB'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_ESC_full.hdf5']='ES_full'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_SRR443884_mESC_2.hdf5']='ES_2'
filenames_sums={}


use_weights=True
split_by_chromosomes=True

def format_plot_axes():
	for i in [1,2]:
		plt.figure(i)
		ax=plt.gca()
		plt.setp(ax.get_xticklabels(), rotation='vertical')
		new_ticks=ax.get_xticks()
		additional_ticks=np.arange(new_ticks[0]+new_ticks[1]/5,new_ticks[1],new_ticks[1]/5)
		new_ticks=additional_ticks.tolist()+new_ticks.tolist()
		ax.set_xticks(new_ticks)
		ax.yaxis.grid(b=True, which='both', color='gray', linestyle='dotted',linewidth=0.5)
		ax.xaxis.grid(b=True, which='both', color='gray', linestyle='dotted',linewidth=0.5)
		if use_weights:
			plt.ylabel('% of contacts')
		else:
			plt.ylabel('Number of contacts')
		plt.xlabel('Distance from current bin')
		plt.title(str(filenames.values())+" weights="+str(use_weights))	
		plt.subplots_adjust(bottom=0.15)
		ax.set_yscale('log',base=2.718281828)
		ax.set_xscale('log',base=2.718281828)

contact_freq={}
contact_freq_sums={}
number_of_bins=0
max_contacts_number = 1
max_contact_length = 0

for i in filenames.keys():  
	print "Reading file "+i+"\n"
	if (i.split('.')[-1]=='hdf5'):
		if (resolution==0): #if we do not know resolution
			raw_heatmap = h5dict.h5dict(i, mode='r') #open heatmap
			resolution = int(raw_heatmap['resolution']) #get the resolution
			del raw_heatmap #close heatmap
		if (genome_db==None): #if we have not initilaize genome before
			genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X']) #do it
		BD = binnedData.binnedData(resolution, genome_db) #now we can initialyze heatmap with defined resolution and genome
		BD.simpleLoad(i, 'HindIII_GM_1')
		number_of_bins=len(BD.dataDict['HindIII_GM_1'])
		
		print "Initiating variables"
		contact_freq[i]={} #For each key "i" containes an array where in position "j" - number of contacts at distance j
		contact_freq_sums[i+'.sum']={} #For each key "i" will contain key "In" with total number of intrachr contacts and Out, with total N of interchr contacts
		for j in range(0,genome_db.chrmCount):
			end0=genome_db.chrmEndsBinCont[j]
			start0=genome_db.chrmStartsBinCont[j]
			contact_freq[i][j]={"In" : [0]*(end0-start0),
						"Out": 0}
			contact_freq_sums[i+'.sum'][j]={"In":[0],
							"Out":0}

			
		print "Calculating contacts"
		for j in range(len(BD.dataDict['HindIII_GM_1'])):
			cur_chr_indx=genome_db.chrmIdxBinCont[j]
			if j % (len(BD.dataDict['HindIII_GM_1']) / 10) == 0:
				print (j / (len(BD.dataDict['HindIII_GM_1']) / 10))*10,"%"
			for k in range(j,len(BD.dataDict['HindIII_GM_1'])):
				if (BD.dataDict['HindIII_GM_1'][j][k] != 0):
					if use_weights:
						contact = BD.dataDict['HindIII_GM_1'][j][k]
					else:
						contact = 1
					if (genome_db.chrmIdxBinCont[k]==cur_chr_indx):
						contact_freq[i][cur_chr_indx]["In"][k-j] += 2*contact
					else:
						contact_freq[i][cur_chr_indx]["Out"] += contact
						contact_freq[i][genome_db.chrmIdxBinCont[k]]["Out"] += contact
					max_contacts_number=max(max_contacts_number,contact)
					max_contact_length=max(max_contact_length,k-j)
		del BD		

contact_length_90=max_contact_length

temps_list=sorted(filenames.keys())

for i in temps_list:
	filenames_sums[i+'.sum']=filenames[i]+'_sum'
	print filenames[i]
	for j in sorted(contact_freq[i].keys()):
		s=sum(contact_freq[i][j]["In"]) # sum of interchromosomal contacts for this chromosome
		print "Chr",j+1,"In=",s,"Out=",contact_freq[i][j]["Out"],"In/out=", \
			float(s)/contact_freq[i][j]["Out"], \
			"In/out(pro MB) ", \
			(float(s)/genome_db.chrmLensBin[j])/(float(contact_freq[i][j]["Out"])/(sum(genome_db.chrmLensBin)-genome_db.chrmLensBin[j]))
		cur_sum=0
		for t in range(1,len(contact_freq[i][j]["In"])):
			cur_sum += contact_freq[i][j]["In"][t]
			if use_weights:
				contact_freq[i][j]["In"][t]=float(contact_freq[i][j]["In"][t])/float(s)
			contact_freq_sums[i+'.sum'][j]["In"].append(float(cur_sum) / s)
							
def contact_freq_total():
	color_numb=0	
	X_values_list=range(1,number_of_bins+1)
	templist={}
	contact_length_90 = number_of_bins
	for i in sorted(filenames.keys()):
		templist[i]=[0]*(number_of_bins+1)
		# Make an average for contact freq for all chromosomes and store it in templist
		for k in range(1,number_of_bins+1):
			s=0
			for j in contact_freq[i].keys():
				if (len(contact_freq[i][j]["In"])>k):
					s+=contact_freq[i][j]["In"][k]
			templist[i][k]=float(s)/len(contact_freq[i].keys())
			if (float(templist[i][k])/templist[i][1]) <= 0.005:
				contact_length_90 = min(contact_length_90,k+1)
		
	
	for i in sorted(filenames.keys()):
			plt.figure(1)
			line_format=colors[color_numb]+markers[color_numb]+'-'
			plt.plot(X_values_list[1:contact_length_90-1],templist[i][1:contact_length_90-1],
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
	del templist

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
	color_numb=0	
	X_values_list=range(1,number_of_bins+1)
	templist={}
	contact_length_90 = max(genome_db.chrmLensBin)
	for i in sorted(filenames_sums.keys()):
		templist[i]=[0]*(number_of_bins+1)
		# Make an average for contact sums for all chromosomes and store it in templist
		for k in range(1,number_of_bins+1):
			s=0
			number_of_contributed_chrs=0
			for j in contact_freq_sums[i].keys():
				if (len(contact_freq_sums[i][j]["In"])>k):
					s+=contact_freq_sums[i][j]["In"][k]
					number_of_contributed_chrs += 1
			if (number_of_contributed_chrs <> 0):
				templist[i][k]=float(s)/number_of_contributed_chrs
			else:
				templist[i][k]=0
#			if (float(templist[i][k])/templist[i][1]) >= 0.999:
#				contact_length_90 = min(contact_length_90,k+1)
		
	
	for i in sorted(filenames_sums.keys()):
			plt.figure(1)
			line_format=colors[color_numb]+markers[color_numb]+'-'
			plt.plot(X_values_list[1:contact_length_90-1],templist[i][1:contact_length_90-1],
					line_format,markersize=1,
					color=colors[color_numb],label=filenames_sums[i])
			color_numb +=1
			if (color_numb >= len(colors)-1 or color_numb>=len(markers)-1):
				color_numb = 0 
	format_plot_axes()
	plt.figure(1)
	plt.legend(loc="center right",bbox_to_anchor=(1.05, 1), borderaxespad=0.)
	filename="Distribution_result_total_"+str(resolution/1000)+"KB_sums"
	if use_weights:
		filename += "_Weighted"
	filename+=".png"
	f = open(filename, "wb")
	plt.savefig(filename,format='png',dpi=300)
	f.close()
	plt.clf()
	del templist


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
contact_freq_by_chromosomes()
contact_freq_total()
contact_sums_by_chromosomes()
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
