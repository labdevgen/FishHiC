print "Starting programm"

import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC
import numpy as np

genome_name='mm9'

genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
domain_res=100000

#domain_coordinates_file_name = "/mnt/storage/home/vsfishman/HiC/data/BASH/Sp_and_Fib_similar_domains_40KB_offset0.txt"
#domain_coordinates_file_name = "/mnt/storage/home/vsfishman/HiC/data/BASH/"#Fib_full2_100KB_all_domains.txt"
bash_input_dir = "/mnt/storage/home/vsfishman/HiC/data/BASH/bash_input/"
annotations_file_path = "/mnt/storage/home/vsfishman/HiC/fasta/mm9/annotation/"
base_folder='/mnt/storage/home/vsfishman/HiC/data/'
#base_filename = 'Fib_full2'

mappability_dic = {}

for base_filename in ['Fib_full2','Sp_full2']:
	domain_coordinates_file_name="/mnt/storage/home/vsfishman/HiC/data/BASH/"+base_filename+"_"+str(domain_res/1000)+"KB_all_domains.txt"
	heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5-raw'
	print "Processing ",domain_coordinates_file_name
	print "Loading file "+heatmap_filepath
	BD = binnedData.binnedData(domain_res, genome_db)
	BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')
	q=BD.dataDict['HindIII_GM_1']
	
	if not genome_db.hasEnzyme():
		genome_db.setEnzyme("HindIII")
	
	chrms = range(genome_db.chrmCount) #numner of chrms in genome
	st = genome_db.chrmStartsBinCont #array of numbers of chrms start (in bins)
	end = genome_db.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)
	
	f=open(domain_coordinates_file_name)
	
	count = 0
	print "Processing domains"
	
	#preparing directory in bash_input
	
	bash_input_dir = bash_input_dir+"/"+domain_coordinates_file_name.split("/")[-1]
	if not os.path.isdir(bash_input_dir):
		os.mkdir(bash_input_dir)
	
	for i in f:
		if len(i.split("\t"))<3:
			print "WARNING: scipping line",i
			continue
		count += 1
		if (count % 100) == 0:
			print count, " processed"
	#		break
		chr=i.strip().split("\t")[0]
		nucleotyde_start=int(i.strip().split("\t")[1])
		nucleotyde_end=int(i.strip().split("\t")[2])
		bin_start = nucleotyde_start / domain_res
		bin_end = nucleotyde_end / domain_res
		with open(bash_input_dir+"/"+domain_coordinates_file_name.split("/")[-1]	+"."+base_filename+"."+chr+"."+str(nucleotyde_start)+"."+str(nucleotyde_end)+".heatmap","w") as f_out:
				#saving a part of heatmap to file
				c = genome_db.label2idx[chr.split('chr')[-1]] #0-based index of current chr
				heatmap=q[st[c]:end[c],st[c]:end[c]][bin_start:bin_end,bin_start:bin_end]
				for x in xrange(len(heatmap)):
					for y in xrange(len(heatmap)):
						if x==y:
							f_out.write("0")
						else:
							f_out.write(str(int(heatmap[x,y])))
						if y != len(heatmap) -1:
							f_out.write("\t")
					if x != len(heatmap)-1:
						f_out.write("\r\n")
		with open(bash_input_dir+"/"+domain_coordinates_file_name.split("/")[-1]+"."+base_filename+"."+chr+"."+str(nucleotyde_start)+"."+str(nucleotyde_end)+".feature","w") as f_out:
			for j in xrange(bin_start,bin_end):
				f_out.write(str(chr)+"\t")
				f_out.write(str(j*domain_res)+"\t")
				f_out.write(str((j+1)*domain_res)+"\t")
				#calculating number or RE sites in locus
				left_enz_index = np.searchsorted(genome_db.rsites[c],j*domain_res,side='left')
				right_enz_index = np.searchsorted(genome_db.rsites[c],(j+1)*domain_res,side='right')
				rsites_number = right_enz_index-left_enz_index
				f_out.write(str(rsites_number)+"\t")
				if rsites_number == 0:
					f_out.write("0.0\t0.0\r\n")
					continue
				#calculating GC content of the locus
				f_out.write(str(genome_db.GCBin[c][j]/100.0)+"\t")
				#calculating mappablility the locus
				if not chr in mappability_dic.keys(): 
					mappability_dic[c]={}
					with open(annotations_file_path+"/"+chr+".an.txt") as f_map_in:
						for line in f_map_in:
							line=line.split("\t")
							try:
								mappability_dic[c][int(line[3])+1]=(float(line[-1])+float(line[-2]))/2.0
							except:
								mappability_dic[c][int(line[3])+1]=0.0
#								print "WARNING: ERROR DURING PARCING LINE",line," (from ",chr,")"
				m=[]
				for t in genome_db.rsites[c][left_enz_index:right_enz_index-1]:
					try:
						m.append(mappability_dic[c][t])
					except:
						print "Chromosome ",chr
						print "rsite=",t
						print "Mappability info not found in ",annotations_file_path+"/"+chr+".an.txt"
						raise "Error!"
				m = np.average(m)
				f_out.write(str(m)+"\r\n")
	f.close()
	print "Done!"
