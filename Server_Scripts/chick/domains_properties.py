import datetime
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
from mirnylib.systemutils import setExceptionHook
setExceptionHook()

from mirnylib import genome
from hiclib import binnedData
from mirnylib import h5dict
import numpy as np
import argparse

def generate_domains_from_file(domains_file):
	domains = np.genfromtxt(domains_file,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
	domains = np.sort(domains,order=["chrm","start"])
	return domains
	
def print_averaged_results_dict(res,printAlldata = False):
	if len(res) == 1:
		for (k,v) in res[0].items():
			output_file.write(str(k)+"\t"+str(v)+"\n")
	else:
		if printAlldata:
			for k in res[0].keys():
				sample = [r[k] for r in res]
				output_file.write(str(k)+"\t"+str(sample)+"\n")
		for k in res[0].keys():
			sample = [r[k] for r in res]
			s = "\t".join(map(str,[k,"average=",np.average(sample),"median=",np.median(sample),"std=",np.std(sample),"min",np.min(sample),"max",np.max(sample)]))
			print s
			output_file.write(s+"\n")


#genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
#				readChrms=[],
#				chrmFileTemplate="N%s.fa")

now = datetime.datetime.now()
print "Starting",str(now)


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")


parser = argparse.ArgumentParser()
parser.add_argument("--hmap")
parser.add_argument("--domains")
parser.add_argument("--genes")
parser.add_argument("--sample",action='store_true')
parser.add_argument("--nolog",action='store_true')
parser.add_argument("--fulllog",action='store_true')
parser.add_argument("--FPKM")
parser.add_argument("--E1")

args = parser.parse_args()

if args.nolog:
	output_file = sys.stdout
else:
	output_file = sys.stdout
	output_file = open("domains_properties_12.txt","a")

output_file.write("------------------------------------------\n")
output_file.write(str(datetime.datetime.now())+"\n")

output_file.write("\t".join(sys.argv)+"\n")

hmap=str(args.hmap).strip()
domains_file=str(args.domains).strip()
if args.sample:
	print "Using shuffle mode to calculate domains statistics"
	bstrap = True
else:
	bstrap = False

print datetime.datetime.now()," loading domains"

if bstrap:
	domains=[generate_domains_from_file(os.path.dirname(domains_file)+"/"+d) for d in os.listdir(os.path.dirname(domains_file)) if os.path.basename(domains_file) in d]
else:
	domains = [generate_domains_from_file(domains_file)]
output_file.write(str(len(domains))+" domain files used\n")

print "Done"


def plot_Contact_drop_depending_on_distance_to_border(domains,distance=1000000,bands_binned= [2,4,6,8,10],colors= ["red","green","blue","black","yellow"]):
	#Contact drop depending on distance to border
	if bstrap:
		print "This function is not designed for bstrap mode"
		print "Skipping function"
		return
	
	print "Contact drop depending on distance to border"
	
	raw_heatmap = h5dict.h5dict(hmap, mode='r') 
	res = int(raw_heatmap['resolution'])
	print "resolution defined by heatmap: ",res

	BD = binnedData.binnedData(res, genome_db)
	print datetime.datetime.now()," loading hmap"
	BD.simpleLoad(hmap, 'heatmap')
	data = BD.dataDict["heatmap"]

	distance_binned = distance / res
	result={}

	for ind,band_binned in enumerate(bands_binned):
		result[band_binned] = {}
		for distance in range(-distance_binned,distance_binned+1):
			result[band_binned][distance] = {}
			result[band_binned][distance]["left"]=[]
			result[band_binned][distance]["right"]=[]
			result[band_binned][distance]["center"]=[]
			for domain in domains:
				chrm = genome_db.label2idx[domain["chrm"]]
				start = int(round(domain["start"] / float(res))) + distance - (band_binned / 2)
				if (start >= 0) and (start + band_binned) < genome_db.chrmLensBin[chrm]:
					start = sum(genome_db.chrmLensBin[0:chrm]) + start
					end = start + band_binned
					result[band_binned][distance]["left"].append(data[start,end])
				
				start = int(round(domain["end"] / float(res))) + distance - (band_binned / 2)
				if (start >= 0) and (start + band_binned) < genome_db.chrmLensBin[chrm]:
					start = sum(genome_db.chrmLensBin[0:chrm]) + start
					end = start + band_binned
					result[band_binned][distance]["right"].append(data[start,end])
				
				domain_length = domain["end"] - domain["start"]
				assert domain_length > 0
				start = int(round((domain["start"] + domain_length/2.) /  float(res))) + distance - (band_binned / 2)
				if (start >= 0) and (start + band_binned) < genome_db.chrmLensBin[chrm]:
					start = sum(genome_db.chrmLensBin[0:chrm]) + start
					end = start + band_binned
					result[band_binned][distance]["center"].append(data[start,end])
				
			result[band_binned][distance]["left"] = np.average(result[band_binned][distance]["left"])
			result[band_binned][distance]["right"] = np.average(result[band_binned][distance]["right"])
			result[band_binned][distance]["center"] = np.average(result[band_binned][distance]["center"])

	print datetime.datetime.now()," Saving pictures"
	for pos in ["left","right","center"]:
		for ind,band_binned in enumerate(bands_binned):
			X = [x*res for x in sorted(result[band_binned].keys())]
			Y = [result[band_binned][x][pos] for x in sorted(result[band_binned].keys())]
			plt.plot(X,Y,label="band="+str(band_binned*res),color=colors[ind],marker="o")
		plt.ylim(ymin=0,ymax=250)
		plt.legend(fontsize="xx-small")
		plt.savefig(hmap+"_"+domains_file.split("/")[-1]+".contact_drop_on_domains_border_"+pos+".png",dpi=300)
		plt.clf()

	print datetime.datetime.now()," Done"


if not bstrap:
	plot_Contact_drop_depending_on_distance_to_border(domains[0],bands_binned=[2])
sys.exit()

#Genes enrichment
print "Reading genes file"
if args.genes is not None:
	genes_file = args.genes.strip()
else:
	#genes_file = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.genes.bed" #CHANGE ME FOR GENES FILE
	#genes_file = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/Gallus_galus.Gallus_gallus-5.0.ncrna.fa.bed" #CHANGE ME FOR ncRNA FILE
	#genes_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/CTCF_Gal5-RepRedv1.bed" #CHANGE ME FOR CTCF FILE
	#genes_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/EBR_input/EBR_chricken100KB_toGalGal5_hglft_genome_7264_fdf1f0.bed" #CHANGE ME FOR EBR FILE
	#genes_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/EBR_input/EBR2_GenomeResearch_hglft_genome_77f7_901430.bed" #CHANGE ME FOR EBR2 FILE
	genes_file = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/Repeats_galGal5.bed" #CHANGE ME FOR TE FILE	


genes = np.genfromtxt(genes_file,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32),('name',"S25")]))
genes = np.sort(genes,order=["chrm","start","end"])
print "Total N of genes:"
print len(genes)

def genes_stistics(domains):
	domains_byChrm = dict([(c,np.atleast_1d(domains[domains["chrm"]==c])) for c in np.unique(domains["chrm"])])
	Ncross, Nwithin, Noutside = 0,0,0
	Nunused_genes = 0

	for chrm in np.unique(genes["chrm"]):
		if not chrm in domains_byChrm.keys():
			if not bstrap:
				print "Warning, no domains on chrm ",chrm," but ",sum(genes["chrm"]==chrm)," genes found there"
			Nunused_genes += sum(genes["chrm"]==chrm)
			continue
		chrm_domains = domains[domains["chrm"]==chrm]
		chrm_genes = genes[genes["chrm"]==chrm]
		stInst_domain_indxs = np.searchsorted(chrm_domains["start"],chrm_genes["start"],side='right')
		stInend_domain_indxs = np.searchsorted(chrm_domains["end"],chrm_genes["start"],side='right')
		endInst_domain_indxs = np.searchsorted(chrm_domains["start"],chrm_genes["end"],side='right')
		endInend_domain_indxs = np.searchsorted(chrm_domains["end"],chrm_genes["end"],side='right')
	
		NotCross =  np.logical_and(np.equal(stInst_domain_indxs,endInst_domain_indxs),np.equal(stInend_domain_indxs,endInend_domain_indxs))
		cross = np.logical_not(NotCross)
		Ncross += sum(cross)

		within = np.logical_and(NotCross,stInst_domain_indxs-endInend_domain_indxs == 1)
		Nwithin += sum(within)
		outside = np.logical_and(NotCross,stInst_domain_indxs-endInend_domain_indxs == 0)
		Noutside += sum(outside)
		
		distances = {"left_in":[],"left_out":[],"right_in":[],"right_out":[]}
		for ind,g in enumerate(chrm_genes):
			if within[ind]:
				d1 = g["start"]-chrm_domains[stInst_domain_indxs[ind]-1]["start"]
				d2 = chrm_domains[stInst_domain_indxs[ind]-1]["end"] - g["end"]
				assert g["start"]>=chrm_domains[stInst_domain_indxs[ind]-1]["start"]
				assert chrm_domains[stInst_domain_indxs[ind]-1]["end"]>=g["end"]
				if d1<d2:
					distances["left_in"].append(d1)
				else:
					distances["right_in"].append(d2)
					
			if outside[ind]:
				if stInst_domain_indxs[ind] != 0 and stInst_domain_indxs[ind] != len(chrm_domains):
					d1=g["start"]-chrm_domains[stInst_domain_indxs[ind]-1]["end"]
					d2=chrm_domains[stInst_domain_indxs[ind]]["start"] - g["end"]
					assert g["start"]>=chrm_domains[stInst_domain_indxs[ind]-1]["end"]
					assert chrm_domains[stInst_domain_indxs[ind]]["start"]>=g["end"]
					if d1<d2:
						distances["right_out"].append(d1)
					else:
						distances["left_out"].append(d2)
				elif stInst_domain_indxs[ind] == 0:
					d = chrm_domains[stInst_domain_indxs[ind]]["start"] - g["end"]
					assert chrm_domains[stInst_domain_indxs[ind]]["start"] >= g["end"]
					distances["left_out"].append(d)
				elif stInst_domain_indxs[ind] == len(chrm_domains):
					d = g["start"]-chrm_domains[stInst_domain_indxs[ind]-1]["end"]
					assert g["start"]>=chrm_domains[stInst_domain_indxs[ind]-1]["end"]
					distances["right_out"].append(d)
				else:
					raise
					
		#DEBUG - slow method of calculation
		# for g in chrm_genes: 
			# o=[overlap_type(g["start"],g["end"],d["start"],d["end"]) for d in chrm_domains]
			# if "cross" in o:
				# Ncross2 += 1
			# elif "within" in o:
				# Nwithin2 += 1
			# elif "outside" in o:
				# assert np.all(np.array(o)=="outside")
				# Noutside2 += 1
			# else:
				# raise
	res = {"cross":Ncross,"within":Nwithin,"outside":Noutside,"Ngenes":len(genes),"Nunused_genes":Nunused_genes}
	for i in distances.keys():
		res[i+"_average_distance_to_gene"] = np.average(distances[i])
	
	if not bstrap:
		print_averaged_results_dict([distances],printAlldata=True)
		
	return res

# print "Calculating genes statistics"
# res = []
# for d in domains:
	# r = genes_stistics(d)
	# res.append(r)
# print_averaged_results_dict(res,printAlldata=True)

def genes_count_depending_on_dist_to_border(domains,
											genes,
											bands=[i*40000 for i in [0.5,1,2]],
											#bands=[i*40000 for i in [1]],											
											max_distance=500000,
#											max_distance=200000,
											step=20000,
#											step=40000,
											colors= ["red","green","blue","black","yellow"],
											FPKMs = None,
											calcFPKM = False,
											calcFPKM_method = np.median,
											file_prefix="",
											StepInDomainPortions = False):
	from intervaltree import Interval, IntervalTree
	chrms_with_no_genes = []
	
	def find_genes_within_interval(chrm,start,end):
		if not chrm in byChrmGenesIntervals:
			if not chrm in chrms_with_no_genes:
				print "[WARNING]: no genes found on chrm ",chrm
				chrms_with_no_genes.append(chrm)
				
			if calcFPKM:
				return []
			else:
				return 0
		if calcFPKM:
			intervalFPKMs = [FPKMs[n.data] for n in byChrmGenesIntervals[chrm][start:end]]
			return intervalFPKMs
#			return max(intervalFPKMs)
		else:
			return len(byChrmGenesIntervals[chrm][start:end])

	byChrmGenes = dict((chrm,genes[genes["chrm"]==chrm]) for chrm in np.unique(genes["chrm"]))
	byChrmGenesIntervals = {}
	print "Generating intervals object"
	for chrm in byChrmGenes.keys():
		 byChrmGenesIntervals[chrm] = IntervalTree( Interval(st,en+1,n) for st,en,n in zip(byChrmGenes[chrm]["start"],byChrmGenes[chrm]["end"],byChrmGenes[chrm]["name"]))
	print "Calculating"
	result = {}

	if StepInDomainPortions == True:
		max_distance = 50
		
		
	for band in bands:
		fulllog_table=np.zeros(shape=(len(domains),len(range(-max_distance,max_distance+step,step))))
		result[band] = { "left":{"x_val":[],"y_val":[],"table":np.zeros_like(fulllog_table)},
						"right":{"x_val":[],"y_val":[],"table":np.zeros_like(fulllog_table)},
						"center":{"x_val":[],"y_val":[],"table":np.zeros_like(fulllog_table)}
						}

		for distance_ind,distance in enumerate(range(-max_distance,max_distance+step,step)):
			result[band]["left"]["x_val"].append(distance)
			result[band]["right"]["x_val"].append(distance)
			result[band]["center"]["x_val"].append(distance)
			if not StepInDomainPortions:
				res_left=[find_genes_within_interval(d["chrm"],d["start"]+distance-band/2.0,d["start"]+distance+band/2.0) for d in domains]
				res_right=[find_genes_within_interval(d["chrm"],d["end"]+distance-band/2.0,d["end"]+distance+band/2.0) for d in domains]
				res_center=[find_genes_within_interval(d["chrm"],(d["end"]+d["start"])/2.+distance-band/2.0,(d["end"]+d["start"])/2.+distance+band/2.0) for d in domains]
			else:
				res_left=[[0]]
				res_right=[[0]]
				res_center=[[0]]
				res_center=[find_genes_within_interval(d["chrm"],
														(d["end"]+d["start"])/2.+distance*(d["end"]-d["start"])/100-band/2.0,
														(d["end"]+d["start"])/2.+distance*(d["end"]-d["start"])/100+band/2.0) for d in domains]
			if calcFPKM:
				result[band]["left"]["table"][:,distance_ind]=[calcFPKM_method(i) if len(i)>0 else 0. for i in res_left]
				result[band]["left"]["y_val"].append(calcFPKM_method(np.concatenate(res_left)))
				
				result[band]["right"]["y_val"].append(calcFPKM_method(np.concatenate(res_right)))
				result[band]["right"]["table"][:,distance_ind]=[calcFPKM_method(i) if len(i)>0 else 0.  for i in res_right]
				
				result[band]["center"]["y_val"].append(calcFPKM_method(np.concatenate(res_center)))
				result[band]["center"]["table"][:,distance_ind]=[calcFPKM_method(i) if len(i)>0 else 0.  for i in res_center]
			else:
				result[band]["left"]["table"][:,distance_ind]=res_left
				result[band]["left"]["y_val"].append(np.average(res_left))
				
				result[band]["right"]["table"][:,distance_ind]=res_right
				result[band]["right"]["y_val"].append(np.average(res_right))
				
				result[band]["center"]["table"][:,distance_ind]=res_center
				result[band]["center"]["y_val"].append(np.average(res_center))
		
			
	
	print datetime.datetime.now()," Saving pictures"

	if not args.nolog:
		output_file.write("genes_count_depending_on_dist_to_border_calcFPKMs="+str(calcFPKM)+"_file_prefix="+file_prefix+"\n")
	
	for pos in result[bands[0]].keys():
		if StepInDomainPortions and pos != "center":
			continue
		if not args.nolog:
			output_file.write("side\t"+pos+"\n")
		for ind,band in enumerate(bands):
			plt.plot(result[band][pos]["x_val"],result[band][pos]["y_val"],label="band="+str(band),color=colors[ind],marker="o")
			if not args.nolog:
				output_file.write("band\t"+str(band)+"\n")
				output_file.write("X:\t"+str(result[band][pos]["x_val"])+"\n")
				output_file.write("Y:\t"+str(result[band][pos]["y_val"])+"\n")
				if args.fulllog:
					for row in xrange(len(result[band][pos]["table"])):
							output_file.write("T:\t"+str(list(result[band][pos]["table"][row]))+"\n")
		plt.legend(fontsize="xx-small")
		plt.savefig(hmap+"_"+domains_file.split("/")[-1]+".genes_drop_on_domains_border_"+pos+"_"+file_prefix+".png",dpi=300)
		plt.clf()
	
	print datetime.datetime.now()," Done!"

def convert_E1_to_genesAndFPKMs_format(E1_file):
	E1array=np.genfromtxt(E1_file,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('value',np.float32)]))
	imitationOfFPKMs={}
	imitationOfGenes=[]
	for ind,val in enumerate(E1array):
		imitationOfGenes.append((val['chrm'],int(val["start"]),int(val["start"])+99999,ind))
		imitationOfFPKMs[str(ind)]=val['value']
	imitationOfGenes = np.sort(np.array(imitationOfGenes,
										dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32),('name','S10')])
										),order=["chrm","start","end"])
	return imitationOfGenes,imitationOfFPKMs
	
def genes_count_crossing_border(domains,genes,max_distance_from_border=5000):
	from intervaltree import Interval, IntervalTree
	byChrmGenes = dict((chrm,genes[genes["chrm"]==chrm]) for chrm in np.unique(genes["chrm"]))
	byChrmGenesIntervals = {}
	for chrm in byChrmGenes.keys():
		byChrmGenesIntervals[chrm] = IntervalTree( Interval(st,en+1,n) for st,en,n in zip(byChrmGenes[chrm]["start"],byChrmGenes[chrm]["end"],byChrmGenes[chrm]["name"]))
	N_cross_left = [len(byChrmGenesIntervals[d["chrm"]][d["start"]-max_distance_from_border:d["start"]+max_distance_from_border]) for d in domains]
	N_cross_right = [len(byChrmGenesIntervals[d["chrm"]][d["end"]-max_distance_from_border:d["end"]+max_distance_from_border]) for d in domains]
	return N_cross_left,N_cross_right

def expression_vs_distance_to_border(domain,genes,FPKMs):
	from intervaltree import Interval, IntervalTree
	byChrmDomains = dict((chrm,domain[domain["chrm"]==chrm]) for chrm in np.unique(domain["chrm"]))
	byChrmDomainsIntervals = {}
	print "Generating intervals object"
	for chrm in byChrmDomains.keys():
		 byChrmDomainsIntervals[chrm] = IntervalTree( Interval(st,en+1,type) for st,en,type in 
								zip(byChrmDomains[chrm]["start"],byChrmDomains[chrm]["end"],["domain"]*len(byChrmDomains[chrm]["end"]))+
								zip(byChrmDomains[chrm]["end"][:-1],byChrmDomains[chrm]["start"][1:],["border"]*(len(byChrmDomains[chrm]["end"])-1))+
								[(0,byChrmDomains[chrm]["start"][0],"border"),(byChrmDomains[chrm]["end"][-1]+1,max(byChrmDomains[chrm]["end"][-1]+2,max(genes["end"])),"border")] )
	
	print "Calculating"
	X = []
	Y = []
	for gene in genes:
		if not gene["chrm"] in byChrmDomainsIntervals.keys():
			continue
		intervals = sorted(byChrmDomainsIntervals[gene["chrm"]][gene["start"]:gene["end"]])
		interval_L,interval_R = intervals[0],intervals[-1]
		assert gene["start"]-interval_L.begin >= 0
		assert interval_R.end-gene["end"] >= 0
		if gene["start"]-interval_L.begin < interval_R.end-gene["end"]:
			if interval_L.data == "domain":
				d = float(gene["start"]-interval_L.begin)
			elif interval_L.data == "border":
				d = -float(gene["start"]-interval_L.begin)
			else:
				raise
		else:
			if interval_R.data == "domain":
				d = float(interval_R.end-gene["end"])
			elif interval_R.data == "border":
				d = -float(interval_R.end-gene["end"])
			else:
				raise
		assert abs(d) <= 50000000
		if abs(d) < 500000:
			X.append(d)
			Y.append(FPKMs[gene["name"]])
		
	plt.clf()
	plt.scatter(X,Y)
	plt.savefig(hmap+"_"+domains_file.split("/")[-1]+".expression_vs_distance_to_border.png")
	plt.clf()
	X2 = range(-1000000,0,20000)+range(0,1000000,20000)
	Y2 = [np.median([y for x,y in zip(X,Y) if left<=x<right]) for left,right in zip(X2[0:-1],X2[1:])]
	X2 = X2[1:]
	plt.plot(X2,Y2)
	plt.savefig(hmap+"_"+domains_file.split("/")[-1]+".expression_vs_distance_to_border2.png")
	plt.clf()
	if not args.nolog:
		output_file.write("expression_vs_distance_to_border2\n")
		output_file.write("X:\t"+str(X2)+"\n")
		output_file.write("Y:\t"+str(Y2)+"\n")

def get_genes_list_at_distance_from_border_or_center(domain,genes,type,distance):
#type = center or border or unsignedBorder
#if border - positive distance for genes inside domain, negative for genes outside
#if unsignedBorder - positive distance to nearest border
	from intervaltree import Interval, IntervalTree
	byChrmDomains = dict((chrm,domain[domain["chrm"]==chrm]) for chrm in np.unique(domain["chrm"]))
	byChrmDomainsIntervals = {}
	print "Generating intervals object"
	for chrm in byChrmDomains.keys():
		 byChrmDomainsIntervals[chrm] = IntervalTree( Interval(st,en+1,type) for st,en,type in 
								zip(byChrmDomains[chrm]["start"],byChrmDomains[chrm]["end"],["domain"]*len(byChrmDomains[chrm]["end"]))+
								zip(byChrmDomains[chrm]["end"][:-1],byChrmDomains[chrm]["start"][1:],["border"]*(len(byChrmDomains[chrm]["end"])-1))+
								[(0,byChrmDomains[chrm]["start"][0],"border"),(byChrmDomains[chrm]["end"][-1]+1,max(byChrmDomains[chrm]["end"][-1]+2,max(genes["end"])),"border")] )
	
	print "Calculating"
	res = []
	for gene in genes:
		if not gene["chrm"] in byChrmDomainsIntervals.keys():
			continue
		intervals = sorted(byChrmDomainsIntervals[gene["chrm"]][gene["start"]:gene["end"]])
		if len(intervals) >= 3:
			continue
		elif type == "center" and len(intervals)==2:
			continue
		elif type == "center":
			center = float(intervals[0].end-intervals[0].begin)/2.
			assert center > 0
			if min(abs(center-float(gene["start"])),abs(float(gene["end"])-center))<=distance:
				res.append(gene["name"])
		else:
			interval_L,interval_R = intervals[0],intervals[-1]
			assert gene["start"]-interval_L.begin >= 0
			assert interval_R.end-gene["end"] >= 0
			if gene["start"]-interval_L.begin < interval_R.end-gene["end"]:
				if interval_L.data == "domain":
					d = float(gene["start"]-interval_L.begin)
				elif interval_L.data == "border":
					d = -float(gene["start"]-interval_L.begin)
				else:
					raise
			else:
				if interval_R.data == "domain":
					d = float(interval_R.end-gene["end"])
				elif interval_R.data == "border":
					d = -float(interval_R.end-gene["end"])
				else:
					raise
			if type == "unsignedBorder":
				if abs(d)<=distance:
					res.append(gene["name"])
			elif type == "border":
				if d<=distance and (d>=0 == distance>=0):
					res.append(gene["name"])
	
	return res
	

def split_domains_to_micro_and_macro(domain):
	from lib_coordinates_converter import *
	converter = Ccoordinates_converter()
	chrmsByType=converter.getChrmsByType()
	chrmsByType2=dict([(chrm,key) for key,val in chrmsByType.iteritems() for chrm in val])
	domainsByType=dict([(key,[]) for key in chrmsByType.keys()])
	from lib_coordinates_converter import *
	converter = Ccoordinates_converter()
	chrmsByType=converter.getChrmsByType()
	chrmsByType2=dict([(chrm,key) for key,val in chrmsByType.iteritems() for chrm in val])
	domainsByType=dict([(key,[]) for key in chrmsByType.keys()])
	for d in domain:
		domainsByType[chrmsByType2[d["chrm"]]].append(d)
	return domainsByType

def compare_micro_and_macro_TADs_sizes(domain):
	output_file.write("compare_micro_and_macro_TADs_sizes\n")
	domainsByType = split_domains_to_micro_and_macro(domain)
	domainsByTypeSizes = dict([(key,[]) for key in domainsByType.keys()])
	for type in domainsByTypeSizes.keys():
		domainsByTypeSizes[type] = [d["end"]-d["start"] for d in domainsByType[type]]
		if not args.nolog:
			hmap_type = "Blood" if "blood" in hmap.lower() else "ChEF"
			domain_type = "Dix" if "dix" in domains_file.lower() else "Arm"
			output_file.write(hmap_type+domain_type+"\t"+type+"\t"+"\t".join(map(str,domainsByTypeSizes[type]))+"\n")
		print type," ",np.average(domainsByTypeSizes[type]),"+/-",np.std(domainsByTypeSizes[type])," ",np.median(domainsByTypeSizes[type])

def genes_count_depending_on_dist_to_border_micro_macro(domains,
											genes,
											file_prefix="",
											**kwargs):
	domainsByType = split_domains_to_micro_and_macro(domains)
	for type in domainsByType.keys():
		genes_count_depending_on_dist_to_border(domainsByType[type],genes,
			file_prefix=file_prefix+"."+type,
			**kwargs)

if args.FPKM is not None: 	#Parse FPKMs file
	print "reading FPKMs file"
	FPKMsArray = np.genfromtxt(args.FPKM,delimiter="\t",dtype=None,skip_header=1)
	gene_names = FPKMsArray["f4"]
	fpkms = FPKMsArray["f9"]
	FPKMs = {}
	for key,val in zip(gene_names,fpkms):
		FPKMs[key] = float(val)
	del FPKMsArray
	count_not_found_genes = 0
	for gene in genes["name"]:
		if not gene in FPKMs.keys():
			count_not_found_genes += 1
			FPKMs[key] = 0.
	print count_not_found_genes," genes not found in FPKM file"
	# genes_count_depending_on_dist_to_border(domains[0],genes,FPKMs=FPKMs,calcFPKM=True,calcFPKM_method=np.median,file_prefix="FPKM_median",bands=[80000])
	genes_count_depending_on_dist_to_border(domains[0],genes,FPKMs=FPKMs,calcFPKM=True,calcFPKM_method=np.median,file_prefix="FPKM_median_portionsStep",bands=[80000],StepInDomainPortions=True,step=10)
if not bstrap:
		print ""
		# genes_count_depending_on_dist_to_border(domains[0],genes,FPKMs=FPKMs,calcFPKM=True,calcFPKM_method=np.average,file_prefix="FPKM_average",bands=[80000])
		
		# print "Defining top and low expressed genes"
		# top = np.percentile(fpkms,75)
		# top_expressed = [key for key,val in FPKMs.iteritems() if val > top]
		# low = np.percentile(fpkms,25)
		# low_expressed = [key for key,val in FPKMs.iteritems() if val < low]
		
		# top_expressed = genes[np.in1d(genes["name"],top_expressed)]
		# low_expressed = genes[np.in1d(genes["name"],low_expressed)]
		
		# print "Done, N_top_expressed = ",len(top_expressed)," N_low_expressed = ",len(low_expressed)
		# genes_count_depending_on_dist_to_border(domains[0],top_expressed,file_prefix="top_expressed",bands=[80000])
		# genes_count_depending_on_dist_to_border(domains[0],low_expressed,file_prefix="low_expressed",bands=[80000])

		# FPKMs2 = {}
		# for gene,fpkm in FPKMs.iteritems():
			# if fpkm >= top:
				# FPKMs2[gene] = 2.
			# elif fpkm < low:
				# FPKMs2[gene] = 0.
			# else:
				# FPKMs2[gene] = 1.
		# genes_count_depending_on_dist_to_border(domains[0],genes,FPKMs=FPKMs2,calcFPKM=True,file_prefix="FPKMsType2Average",calcFPKM_method=np.average,bands=[80000])
#		expression_vs_distance_to_border(domains[0],genes,FPKMs)
	# print datetime.datetime.now()," Calculating genes drop on domains border"
	# genes_count_depending_on_dist_to_border(domains[0],genes)
	# print datetime.datetime.now()," Done"

	# print datetime.datetime.now()," Calculating genes crossing domains border"
	# N_cross_left,N_cross_right = genes_count_crossing_border(domains[0],genes)

	# if args.nolog:
		# print "Left average: ",np.average(N_cross_left)
		# print "Right average: ",np.average(N_cross_right)
	# else:
		# output_file.write("genes_count_crossing_border\n")
		# output_file.write("Left:\t"+str(N_cross_left)+"\n") 
		# output_file.write("Right:\t"+str(N_cross_right)+"\n")
	# print datetime.datetime.now()," Done"
	
	# print datetime.datetime.now()," Generating list of genes for GO"
	# border_genes=get_genes_list_at_distance_from_border_or_center(domains[0],genes,"unsignedBorder",10000)
	# center_genes=get_genes_list_at_distance_from_border_or_center(domains[0],genes,"center",10000)
	# output_file2 = open("domain_properties.GO.txt","a")
	# output_file2.write("------------------------------------------\n")
	# output_file2.write(str(datetime.datetime.now())+"\n")
	# output_file2.write("\t".join(sys.argv)+"\n")
	# output_file2.write("get_genes_list_at_distance_from_border_or_center: border\n")
	# output_file2.write("\n".join(border_genes)+"\n")
	# output_file2.write("get_genes_list_at_distance_from_border_or_center: center\n")
	# output_file2.write("\n".join(center_genes)+"\n")
	# output_file2.write("Done\n")
	# output_file2.close()


genes_count_depending_on_dist_to_border(domains[0],genes,FPKMs=None,calcFPKM=False,file_prefix="TE_elements",bands=[80000],step=40000,max_distance=500000)

#output_file.write(str(datetime.datetime.now())+" Done!\n")
#compare_micro_and_macro_TADs_sizes(domains[0])
#genes_count_depending_on_dist_to_border(domains[0],genes,file_prefix="E1",bands=[20000,40000,80000,120000],FPKMs=None,calcFPKM=False)
#genes_count_depending_on_dist_to_border_micro_macro(domains[0],genes,file_prefix="ncRNAgenes",bands=[20000,40000,80000,120000],FPKMs=None,calcFPKM=False)

if args.E1 is not None:
	imitationOfGenes,imitationOfFPKMs = convert_E1_to_genesAndFPKMs_format(args.E1.strip())
	genes_count_depending_on_dist_to_border(domains[0],imitationOfGenes,file_prefix="E1",bands=[100000],step=50000,FPKMs=imitationOfFPKMs,calcFPKM=True,calcFPKM_method=np.average)

output_file.close()