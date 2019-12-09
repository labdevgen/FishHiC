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

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--hmap",help='path to heatmap file')
parser.add_argument("--reg",help='regions\ne.g. chr1:1-100,chr2:1-10000\nor"all" to save by chromosome pics')
parser.add_argument("--level",default="chr",help='"chr" for chromosome level, "contig" for congig level')
parser.add_argument("--domains",default="",help='domains file to use for annotation,coud be several files separated with ","')
parser.add_argument("--colors",default="",help='colors for domains,should be same number of colors as number of domains')
parser.add_argument("--out",default="heatmap_plots",help='output dir for plots')

args = parser.parse_args()
if args.level == "contig":
	genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")
elif args.level == "chr":
	genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")
	#genome_db = Genome('/mnt/storage/home/vsfishman/HiC/fasta/mm9', readChrms=['#', 'X'])

hm_file = args.hmap
figure_path = args.out+"/"+hm_file.split("/")[-1]

if args.domains != "":
	domains=args.domains.split(",")
	colors=args.colors.split(",")
else:
	domains=[]
	colors=[]

assert len(domains)==len(colors)
domains=[np.genfromtxt(domain_file,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
					for domain_file in domains]
resolution = extractResolutionFromFileName(hm_file)
print "Using resolution ",resolution

data_dict = {} # a structure to keep chrms arrays
data_dict_second_hmap = {} # a structure to keep chrms arrays
				
for chrm,array in enumerate(get_chromosomes(hm_file,genome_db,resolution)):
	data_dict[chrm] = array

if args.reg == "all":
	args.reg = ",".join([i+":0-"+str(genome_db.chrmLens[genome_db.label2idx[i]]) for i in genome_db.chrmLabels 
												if genome_db.chrmLens[genome_db.label2idx[i]] > resolution*10])
	
for i in args.reg.split(","):
		chrm = i.split(":")[0]
		nt_start,nt_end = map(int,i.split(":")[1].split("-"))
		title = "_".join(map(str,[chrm,nt_start,nt_end]))
		chrmLabel = chrm
		chrm = genome_db.label2idx[chrm]
		data = data_dict[chrm]
		left_border = nt_start/resolution
		right_border = nt_end/resolution
		to_plot = np.log(data[left_border:right_border,left_border:right_border])
			

		# plt.plot([distance/resolution,right_border-left_border-(distance)/resolution,(right_border-left_border)-distance/resolution],
					# [distance/resolution,distance/resolution,-(distance/resolution)+right_border-left_border],
					# ls="dashed",color='k')
			
		plt_inches = plt.figure().get_size_inches() # number of points in plot, array [X,Y]
		inches_per_dot = float(min(plt_inches)) / len(to_plot) #number of data-points / inch
		points_per_dot = max(0.1,(inches_per_dot*72)/4) #number of data-points per "point" (1 inch = 72 points in matplotlib)
		print chrm,points_per_dot
		#find_domains_inside_region
		for domain,domain_color in zip(domains,colors):
			if ("chr" in domain[0]["chrm"]) and (not "CHR" in chrmLabel.upper()):
				addition = "CHR"
			else:
				addition = ""
			domains_in_region = [d for d in domain if d["start"]>=nt_start and 
														d["end"]<=nt_end and 
														d["chrm"].upper()==addition+chrmLabel.upper()]
			for d in domains_in_region:
					moving_constant = nt_start
					st = (d["start"] - moving_constant)/resolution
					end =(d["end"] - moving_constant)/resolution
					plt.plot([st,st,end],
							[st,end,end],
							ls="solid",linewidth=points_per_dot,color=domain_color)
		plotting.plot_matrix(to_plot,cmap='OrRd')
		if right_border-left_border > 10:
			tick_coeff = 10
		else:
			tick_coeff = 1
		plt.yticks(list(range(0,right_border-left_border,((right_border-left_border)/tick_coeff))))
		plt.xticks(list(range(0,right_border-left_border,((right_border-left_border)/tick_coeff))))
		plt.gca().set_yticklabels(["%.5g" % i for i in range(left_border*resolution,
															right_border*resolution,
															((right_border-left_border)/tick_coeff)*resolution)])
		plt.gca().set_xticklabels(["%.5g" % i for i in range(left_border*resolution,
															right_border*resolution,
															((right_border-left_border)/tick_coeff)*resolution)],
								rotation="vertical")

		plt.title(title)
		plt.subplots_adjust(bottom=0.15)
		current_figure_path = figure_path+"."+title+".png"
		f = open(current_figure_path, "wb")
		plt.savefig(current_figure_path,dpi=900)
		f.close()
		plt.clf()