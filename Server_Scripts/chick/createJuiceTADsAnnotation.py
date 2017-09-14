import sys

try:
		import numpy as np
except:
		sys.path = ["/usr/lib64/python2.7/site-packages"]+sys.path
		import numpy as np


import os
import gzip
from mirnylib.systemutils import setExceptionHook
setExceptionHook()

if len(sys.argv) < 4:
	print "-------------Usage:----------"
	print sys.argv[0]," path_do_directory_with_domains mode color (resolution)"
	print "mode: a = armatus domains, d = dixon domains t = TADtree domains"
	print "for TADtree mode please also provide resolution"
	print "color: RGB color, e.g. 255,0,0"
	sys.exit()


dir=sys.argv[1].strip()
mode = sys.argv[2]
color = sys.argv[3]
if dir.endswith("/"):
	dir = dir[:-1]

out = open(dir+"/"+dir.split("/")[-1]+".jucebox_domains.annotation","w")
print "Saving file",dir+"/"+dir.split("/")[-1]+".jucebox_domains.annotation"

if mode == "a":
	end = ".gz.domains.consensus.txt"
elif mode == "d":
	from mirnylib import genome
	genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")
	filter = lambda x: x.endswith(".final_domains")
elif mode == "t":
	resolution = int(sys.argv[4])
	import numpy as np
	import numpy.lib.recfunctions
	from mirnylib import genome
	from intervaltree import Interval, IntervalTree
	genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")
	filter = lambda x: x.startswith("chr")
	def process_TADtree(fname):
		dup=np.atleast_1d(np.genfromtxt(dir+"/"+fname+"/proportion_duplicates.txt",dtype=None,skip_header=1))
		N_index = min(len(dup)-1,np.searchsorted(dup["f1"],0.02))
		N = dup["f0"][N_index]
		print "For chr",fname," using N=",N,"; total N=",dup["f0"][-1]
		domains = np.atleast_1d(np.genfromtxt(dir+"/"+fname+"/"+N+".txt",dtype=None,skip_header=1))
		length = (1./(domains["f2"]-domains["f1"]))
		domains = numpy.lib.recfunctions.append_fields(domains,["l"],[length],usemask=False)
		domains = np.sort(domains,order=["f1","l"])
		domains["f1"] = (domains["f1"]-1)*resolution# - resolution / 2
		domains["f2"] = (domains["f2"]-1)*resolution# - resolution / 2
		assert np.all(domains["f1"]>=0)
		assert np.all(domains["f2"]>0)
		levels = []
		domainsTree = IntervalTree()
		for i in domains:
			st = i["f1"]
			end = i["f2"]
			levels.append(len(domainsTree[st+1:end-1]))
			domainsTree[st+1:end-1] = 1
		assert len(levels)==len(domains)
		return levels,domains

for fname in [f for f in os.listdir(dir+"/") if filter(f)]:
	print "Processing ",fname
	if mode=="t":
		levels,domains = process_TADtree(fname)
		print "Max Level = ",max(levels)
		for i,level in zip(domains,levels):
			out.write("\t".join(map(str,[i["f0"],i["f1"],i["f2"],\
						i["f0"],i["f1"],i["f2"],\
						color,"Level_"+str(level)+"_"+dir.split("/")[-1]]))+"\n")
		continue
	file = open(dir+"/"+fname)
	for line in file:
		if "Nchr" in line:
			line = line.replace("Nchr","chr",1)
		if mode == "d":
			line = line.strip().split()
			line[0] = genome_db.idx2label[int(line[0].split("chr")[-1])]
			assert len(line) == 3
			line = "\t".join(line)
		out.write(line.strip()+"\t"+line.strip()+"\t"+color+"\t"+dir.split("/")[-1]+"\n")
	file.close()
	
out.close()
print "Done!"
