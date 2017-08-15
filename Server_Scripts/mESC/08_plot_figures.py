import sys
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
from  getIntraChrHeatmaps import get_chromosomes,extractResolutionFromFileName
import subprocess

def executeBashCommand(bashCommand,doPrint=True):
  print "executing command %s \n" % (bashCommand)
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  if doPrint:
    print p.communicate()[0]
    return ""
  else:
    return p.communicate()[0]


genomeName = "mm10"
genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/"+genomeName,readChrms=["#","X","Y"])

chainFail = "/mnt/storage/home/vsfishman/HiC/fasta/mm10/mm10ToMm9.over.chain"
#hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res25k_hiRes.hm.IC.cis.hdf5"
#hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res100k_bychr.hm.IC.hdf5"
hm_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5"

#results_file = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res25k_hiRes.hm.IC.cis.hdf5.good_examples.5.6.7.txt"
results_files = ["/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/"+fname for fname in os.listdir("/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/")
					if os.path.isfile("/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/"+fname) and ("good_examples" in fname) and (not ("criteias") in fname) and (hm_file.split("/")[-1] in fname)] 
#results_files = ["/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5.good_examples.4.6.txt"]
print "Processing files:\n"+"\n".join(results_files)
threshold=99
distance = 100			

#points = generate_points_from_file(hm_file+".good_examples.txt",threshold=95)
#points = generate_points_from_file("/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5.good_examples4.5.6.7.txt",
#									threshold=95)

def generate_points_from_file(fname,threshold=0.):
	points = []
	f=open(fname,"r")
	f.readline() #skip a header row
	for line in f:
		line=line.strip().split()
		if len(line)<=3:
			print "point not found in line\n","\t".join(line)
		if float(line[2])>=threshold:
			chrm = int(line[0])
			nt = (float(line[1]) - float(resolution)/2.)
			assert (nt % resolution) == 0
			points.append((chrm,nt))
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

assert len(all_result_files_points) > 0
	
data_dict = {} # a structure to keep chrms arrays
if len(np.unique([p[0] for p in points for points in all_result_files_points.values()]))==1:
	data_dict[points[0][0]] = get_chromosomes(hm_file,genome_db,resolution,chrNumb=points[0][0])
else:
	for chrm,array in enumerate(get_chromosomes(hm_file,genome_db,resolution)):
		data_dict[chrm] = array
	
for results_file in all_result_files_points:
	points = all_result_files_points[results_file]
	figs_dir = results_file+".tr_"+str(threshold)+".figs/"
	for i in points:
		chrm,nt = i
		bin = int(nt/resolution)
		data = data_dict[chrm]
		#data is a matrix (2D numpu array)
		left_border = max(0,bin-distance)
		right_border = min(len(data),bin+distance)
		to_plot = np.log(data[left_border:right_border,left_border:right_border])
		temp = np.array(to_plot).flatten()
		temp = temp[temp.nonzero()]
		to_plot[bin-left_border,bin-left_border]=np.min(temp[np.isfinite(temp)])
		figure_path = figs_dir+str(resolution)+"_"+str(chrm)+"_"+str(bin)+".png"
		print "saving ",figure_path
		plotting.plot_matrix(to_plot)
	
		mm9pos = convert_points(chrm,left_border*resolution,right_border*resolution,"/mnt/storage/home/vsfishman/HiC/fasta/mm10/mm10ToMm9.over.chain")
		if mm9pos != "N/A":
			exect_mm9pos = (int(mm9pos.split(":")[-1].split("-")[0])+int(mm9pos.split(":")[-1].split("-")[1]))/2
		else:
			exect_mm9pos = "N/A"
		plt.title("mm10 chr"+str(chrm+1)+":"+str(left_border*resolution)+"-"+str(right_border*resolution) + 
											" ("+str(bin*resolution) + ")\nmm9:"+mm9pos+"("+str(exect_mm9pos)+")")
		plt.subplots_adjust(bottom=0.15)	
		f = open(figure_path, "wb")
		plt.savefig(figure_path,dpi=300)
		f.close()
		plt.clf()