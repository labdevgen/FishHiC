import sys

try:
		import numpy as np
except:
		sys.path = ["/usr/lib64/python2.7/site-packages"]+sys.path
		import numpy as np


import os
import gzip

if len(sys.argv) != 4:
	print "-------------Usage:----------"
	print sys.argv[0]," path_do_directory_with_domains mode color"
	print "mode: a = armatus domains, d = dixon domains"
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
	end = ".final_domains"

for fname in [f for f in os.listdir(dir+"/") if f.endswith(end)]:
#	print "Processing ",fname
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