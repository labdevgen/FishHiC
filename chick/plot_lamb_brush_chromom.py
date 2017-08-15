import sys
import math
try: 
	import numpy as np 
except:
	sys.path = ["/usr/lib64/python2.7/site-packages"]+sys.path
	import numpy as np
print "Numpy inported!"
from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome
from hiclib import highResBinnedData
from mirnylib import h5dict
from mirnylib import plotting
from lib_coordinates_converter import *

converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
converter.prepare_acc2chrmDict()							
converter.create_agp_dict()							


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")
from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")


genome_db_contig = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")
genome_db_chrmLevel = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")
genome_db_contig = genome_db_chrmLevel

#hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII-100k.hm.IC"
#second_hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII-100k.hm.IC"

########################WRITE YOUR HEATMAP HERE########################
hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.IC"
domains_files_Arm = "mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.gzipped_matrix/ChEF-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation"
domains_files_Dix = "/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation"

second_hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.IC"
second_domains_files_Arm = "mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.gzipped_matrix/Blood-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation"
second_domains_files_Dix = "/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation"


resolution = extractResolutionFromFileName(hm_file)
assert resolution == extractResolutionFromFileName(second_hm_file)

data_dict = {} # a structure to keep chrms arrays
data_dict_second_hmap = {} # a structure to keep chrms arrays
				
for chrm,array in enumerate(get_chromosomes(hm_file,genome_db_contig,resolution)):
	data_dict[chrm] = array
for chrm,array in enumerate(get_chromosomes(second_hm_file,genome_db_contig,resolution)):
	data_dict_second_hmap[chrm] = array

#chrm,nt_start,nt_end,title

#######################################ADD REGIONS TO PLOT HERE
###EXAMPLE:
###points = [
#["chr4",40100000,44000000,"Chromomer_Anja"],
#["chr5",4000,440000,"Chromomer_Anja"]
#]

points = [["chr15",7200000,7280000,"Gamatu_Gene"]]
#points = [["chr15",6200000,8280000,"Gamatu_Gene"]]
distance = 1000000

for i in points:
	
		chrm,nt_start,nt_end,title = i
		#position = converter.chrmTOcontig(chrm,nt_start)
		#chrm1,nt_start = position.chrm,position.nt
		#position = converter.chrmTOcontig(chrm,nt_start)
		#chrm2,nt_end = position.chrm,position.nt
		#assert chrm1==chrm2
		#chrm=chrm1
		
		chrmLabel = chrm
		chrm = genome_db_contig.label2idx[chrm]
		
		for data,hmpath,domains in zip([data_dict[chrm],data_dict_second_hmap[chrm]],
										[hm_file,second_hm_file],
										[[domains_files_Arm,domains_files_Dix],
										[second_domains_files_Arm,second_domains_files_Dix]]):
			left_border = (nt_start-distance)/resolution
			right_border = (nt_end+distance)/resolution
			to_plot = np.log(data[left_border:right_border,left_border:right_border])
			
			plotting.plot_matrix(to_plot)
			
			plt.plot([distance/resolution,right_border-left_border-(distance)/resolution,(right_border-left_border)-distance/resolution],
					[distance/resolution,distance/resolution,-(distance/resolution)+right_border-left_border],
					ls="dashed",color='k')
			
			#find_domains_inside_region
			for domain_file,domains_color in zip(domains,["black","white"]):
				domain = np.genfromtxt(domain_file,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
				domain = np.sort(domain,order=["chrm","start"])
				domains_in_region = [d for d in domain if d["start"]>=(nt_start-distance) and d["end"]<=(nt_end+distance) and d["chrm"].upper()==chrmLabel.upper()]
				for d in domains_in_region:
					moving_constant = nt_start-distance
					st = (d["start"] - moving_constant)/resolution
					end =(d["end"] - moving_constant)/resolution
					plt.plot([st,st,end],
							[st,end,end],
							ls="solid",color=domains_color)
			plt.title(title)
			plt.subplots_adjust(bottom=0.15)
			figure_path = hmpath+"."+title+".png"
			f = open(figure_path, "wb")
			plt.savefig(figure_path,dpi=300)
			f.close()
			plt.clf()