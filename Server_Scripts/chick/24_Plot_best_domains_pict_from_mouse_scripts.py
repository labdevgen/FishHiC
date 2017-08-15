print "Starting programm"
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")
sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/")
from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName
from lib_matrix_operations import *

from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName

import scipy.stats
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC
import math

genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

hm_file = sys.argv[1]
figure_path = "/mnt/storage/home/vsfishman/HiC/pics/"

domain_res = extractResolutionFromFileName(hm_file)
all_chrms=get_chromosomes(hm_file,genome_db,domain_res)

chrms = range(genome_db.chrmCount) #numner of chrms in genome
st = genome_db.chrmStartsBinCont #array of numbers of chrms start (in bins)
end = genome_db.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)


def plot_one_chr_fragment(matrix,figure_path,chr_name,ch_start,ch_end): 
	print "Plotting picture"
	ch_start = ch_start / domain_res
	ch_end = ch_end / domain_res
	# domain_st = domain_st / domain_res
	# domain_end = domain_end / domain_res
	i=genome_db.label2idx[chr_name]
	q2_0=matrix[i]
	q2_0=q2_0[ch_start:ch_end,ch_start:ch_end]
#	np.savetxt("test_for_plot.txt",q2_0)
	plotting.plot_rotated_matrix(np.log(q2_0),0,cmap='OrRd')
	print ch_start,ch_end#,domain_st,domain_end
#	print domain_st-ch_start,domain_end-ch_start
#	plt.axvline(x=domain_st-ch_start)
#	plt.axvline(x=domain_end-ch_start)
	plt.subplots_adjust(bottom=0.15)
	fp1=figure_path+"_"+chr_name+"_bin_"+str(ch_start)+"_to_"+str(ch_end)+".rotated.png"
	print "Saving figure "+fp1
	f = open(fp1, "wb")
	plt.savefig(fp1,dpi=900)
	f.close()
	
	plt.clf()
	plotting.plot_matrix(np.log(q2_0),cmap='OrRd')
	fp1=figure_path+"_"+chr_name+"_bin_"+str(ch_start)+"_to_"+str(ch_end)+".png"
	print "Saving figure "+fp1
	f = open(fp1, "wb")
	plt.savefig(fp1,dpi=900)
	f.close()


#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr2",9200000,9600000)
plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr2",31500000,34000000)
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr2",7600001,9600000)
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr4",0,genome_db.chrmLens[genome_db.label2idx["chr4"]])
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr12",0,genome_db.chrmLens[genome_db.label2idx["chr12"]])
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr9",1000000,9000000)
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr1",16250000,24250000)
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr6",11400000,20400000)
#plot_one_chr_fragment(all_chrms,figure_path+hm_file.split("/")[-1],"chr4",46600000,55600000)