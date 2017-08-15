#SCRIPT USED TO FIX A BUG produced by convertGenomes.py
#Bug is already fixed in code of convertGenomes.py, so basicly this script should not be used any more


import os
import numpy as np
from mirnylib.h5dict import h5dict

from mirnylib.systemutils import setExceptionHook
setExceptionHook()
from mirnylib import genome

genome_db_chrmLevel = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

new_genome_dict={"genome":{"label2idx":genome_db_chrmLevel.label2idx,"idx2label":genome_db_chrmLevel.idx2label}}

dir_to_convert = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/"
files_to_convert = []

for dirpath,dirnames,filenames in os.walk(dir_to_convert):
	if "completed" in filenames:
		files_to_convert += [os.path.join(dirpath,f) for f in filenames if f.endswith(".updGenome")]
	elif sum([f.endswith(".updGenome") for f in filenames])>0:
		print "Ommiting files in directory ",dirpath
print "Converting files:\n","\n".join(files_to_convert)

for fname in files_to_convert:
	data = h5dict(fname, mode="r+")
	data["misc"]=new_genome_dict
	data.update()
	print fname," updated"