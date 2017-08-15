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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import subprocess

def executeBashCommand(bashCommand,doPrint=True):
  print "executing command %s \n" % (bashCommand)
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  if doPrint:
    print p.communicate()[0]
    return ""
  else:
    return p.communicate()[0]


sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/mESC")
from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName

genome_db_contig = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")

hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII-100k.hm.IC"
second_hm_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII-100k.hm.IC"

#results_files = ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII-100k.hm.IC.good_examples.10.12.14.txt"]
results_files = ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/"+fname 
					for fname in os.listdir("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/")
					if os.path.isfile("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/"+fname) and 
					("good_examples" in fname) and (not ("criteias") in fname) and (hm_file.split("/")[-1] in fname)] 
#results_files = ["/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5.good_examples.4.6.txt"]
print "Processing files:\n"+"\n".join(results_files)
threshold=95
distance = 100

#points = generate_points_from_file(hm_file+".good_examples.txt",threshold=95)
#points = generate_points_from_file("/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5.good_examples4.5.6.7.txt",
#									threshold=95)

def generate_points_from_file(fname,threshold=0.):
	points = []
	f=open(fname,"r")
	header=f.readline().strip().split() #skip a header row
	for line in f:
		line=line.strip().split()
		if len(line)<=5:
			print "point not found in line\n","\t".join(line)
		if float(line[4])>=threshold:
			chrm = genome_db_contig.label2idx[line[0][1:]]
			nt = (float(line[1]) - float(resolution)/2.)
			assert (nt % resolution) == 0
			best_dist = map(abs,map(float,line[5:]))
			best_dist = best_dist.index(max(best_dist))
			best_dist = int(header[5+best_dist])
			assert distance > best_dist + 1
			points.append((chrm,nt,best_dist,float(line[4]),line[2],float(line[3])))
	f.close()
	return points

def convert_points(chr,start,end,chainFail):
	temp_file = open("tempOldCoorfinatesFile.txt","w")
	if chr == 19:
			chrm = "chrX"
	elif chr == 20:
			chrm = "chrY"
	else:
			chrm = "chr"+str(chr+1)
	temp_file.write(chrm+":"+str(int(start))+"-"+str(int(end)))
	temp_file.close()
	executeBashCommand("liftOver -positions tempOldCoorfinatesFile.txt "+chainFail+" tempOldCoorfinatesFile.conv.txt tempOldCoorfinatesFile.umap.txt")
	try:
		f=open("tempOldCoorfinatesFile.umap.txt")
		for line in f:
			if line[0]=="#":
				continue
			else:
				f.close()
				return "N/A"
	except:
		pass
	f=open("tempOldCoorfinatesFile.conv.txt")
	for line in f:
		if line[0]=="#":
			continue
		else:
			f.close()
			return line.strip()

all_result_files_points = {}
for results_file in results_files:
	resolution = extractResolutionFromFileName(hm_file)
	points = generate_points_from_file(results_file,
									threshold=threshold)
	if len(points)==0:
		print "No points relevant for threshold found in a file",results_file
		continue
	all_result_files_points[results_file] = points
	figs_dir = results_file+".tr_"+str(threshold)+".figs/"
	if os.path.exists(figs_dir):
		import shutil
		shutil.rmtree(figs_dir)

	os.mkdir(figs_dir)
	
final_figs_dir = hm_file+".FISH.figs/"
if os.path.exists(final_figs_dir):
		import shutil
		shutil.rmtree(final_figs_dir)

os.mkdir(final_figs_dir)
	

assert len(all_result_files_points) > 0
	
data_dict = {} # a structure to keep chrms arrays
data_dict_second_hmap = {} # a structure to keep chrms arrays
if len(np.unique([p[0] for p in points for points in all_result_files_points.values()]))==1:
	data_dict[points[0][0]] = get_chromosomes(hm_file,genome_db_contig,resolution,chrNumb=points[0][0])
	data_dict_second_hmap[points[0][0]] = get_chromosomes(second_hm_file,genome_db_contig,resolution,chrNumb=points[0][0])
else:
	for chrm,array in enumerate(get_chromosomes(hm_file,genome_db_contig,resolution)):
		data_dict[chrm] = array
	for chrm,array in enumerate(get_chromosomes(second_hm_file,genome_db_contig,resolution)):
		data_dict_second_hmap[chrm] = array
		
final_results = open(hm_file+".FISH_regions.txt","w")
for results_file in all_result_files_points:
	points = all_result_files_points[results_file]
	figs_dir = results_file+".tr_"+str(threshold)+".figs/"
	for i in points:
		chrm,nt,best_dist,point_treashold,chrmChrmLevel,ntChrmLevel = i
		bin = int(nt/resolution)
		data = data_dict[chrm]
		
		left = data_dict_second_hmap[chrm][bin,bin-best_dist]
		right = data_dict_second_hmap[chrm][bin,bin+best_dist]
		ratio = data[bin,bin-best_dist] / data[bin,bin+best_dist]

		if (point_treashold >= 99) or (left*right != 0 and math.log(left/right)*math.log(ratio) < 0):
			#data is a matrix (2D numpu array)
			avaliable_distance = min(distance,bin,len(data)-bin-1)
			assert avaliable_distance > 0
			assert best_dist <= avaliable_distance
			left_border = bin-avaliable_distance
			right_border = bin+avaliable_distance+1
			bin_realative_number = avaliable_distance
			assert bin_realative_number>0
			to_plot = np.log(data[left_border:right_border,left_border:right_border])
			temp = np.array(to_plot).flatten()
			temp = temp[temp.nonzero()]
			to_plot[bin-left_border,bin-left_border]=np.min(temp[np.isfinite(temp)])

			plotting.plot_matrix(to_plot)
			plt.scatter([bin_realative_number-best_dist,bin_realative_number],
						[bin_realative_number,bin_realative_number+best_dist],
						facecolors='none', edgecolors='r')
			#Draw triangles to visualize regions
			plt.plot([bin_realative_number-best_dist,bin_realative_number,bin_realative_number,bin_realative_number+best_dist,bin_realative_number+best_dist],
					[bin_realative_number-best_dist,bin_realative_number-best_dist,bin_realative_number,bin_realative_number,bin_realative_number+best_dist],
					ls="dashed") 

			plt.title("N"+genome_db_contig.idx2label[chrm] + " " + str(nt+float(resolution)/2.)+"+/-"+str(best_dist*100)+"kb")
			plt.subplots_adjust(bottom=0.15)
		
			figure_path = final_figs_dir+str(resolution)+"_"+str(chrm)+"_"+str(bin)+".png"
			f = open(figure_path, "wb")
			plt.savefig(figure_path,dpi=300)
			f.close()
			plt.clf()
		else: 
			figure_path = figs_dir+str(resolution)+"_"+str(chrm)+"_"+str(bin)+".png"

		if left*right != 0 and math.log(left/right)*math.log(ratio) < 0:
			print "Suitable!"
			data2 = data_dict_second_hmap[chrm]
			to_plot2 = np.log(data2[left_border:right_border,left_border:right_border])
			temp = np.array(to_plot2).flatten()
			temp = temp[temp.nonzero()]
			to_plot2[bin-left_border,bin-left_border]=np.min(temp[np.isfinite(temp)])
			to_plot2[bin-left_border]
			figure_path = final_figs_dir+str(resolution)+"_"+str(chrm)+"_"+str(bin)+"_second_sample.png"
			print "saving ",figure_path
			plotting.plot_matrix(to_plot2)
			plt.scatter([bin_realative_number-best_dist,bin_realative_number],
						[bin_realative_number,bin_realative_number+best_dist],
						facecolors='none', edgecolors='r')

			#Draw triangles to visualize regions
			plt.plot([bin_realative_number-best_dist,bin_realative_number,bin_realative_number,bin_realative_number+best_dist,bin_realative_number+best_dist],
					[bin_realative_number-best_dist,bin_realative_number-best_dist,bin_realative_number,bin_realative_number,bin_realative_number+best_dist],
					ls="dashed") 
			
			plt.title("N"+genome_db_contig.idx2label[chrm] + " " + str(nt+float(resolution)/2.)+"+/-"+str(best_dist*100)+"kb")
			plt.subplots_adjust(bottom=0.15)
			f = open(figure_path, "wb")
			plt.savefig(figure_path,dpi=300)
			f.close()
			plt.clf()
			final_results.write("\t".join(map(str,["N"+genome_db_contig.idx2label[chrm],
													best_dist*resolution,
													nt,
													nt-best_dist*resolution,
													nt+best_dist*resolution,
													chrmChrmLevel,
													ntChrmLevel,
													ratio,
													"Both",
													left/right,
													str(resolution)+"_"+str(chrm)+"_"+str(bin)+".png"]))+"\n")
		else:
			if point_treashold >= 99:
				final_results.write("\t".join(map(str,["N"+genome_db_contig.idx2label[chrm],
													best_dist*resolution,
													nt,
													nt-best_dist*resolution,
													nt+best_dist*resolution,
													chrmChrmLevel,
													ntChrmLevel,
													ratio,
													hm_file.split("/")[-1].split("-")[0],
													"n/a",
													str(resolution)+"_"+str(chrm)+"_"+str(bin)+".png"]))+"\n")

final_results.close()