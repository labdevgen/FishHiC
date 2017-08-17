import os
import numpy as np
from mirnylib.h5dict import h5dict
from lib_coordinates_converter import *

from mirnylib.systemutils import setExceptionHook
setExceptionHook()
from mirnylib import genome
genome_db_contigLevel = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")
genome_db_chrmLevel = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")


dir_to_convert = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/"
files_to_convert = []

for dirpath,dirnames,filenames in os.walk(dir_to_convert):
	if "completed" in filenames:
		files_to_convert += [os.path.join(dirpath,f) for f in filenames if f.endswith(".hdf5")]
	elif sum([f.endswith(".hdf5") for f in filenames])>0:
		print "Ommiting files in directory ",dirpath
print "Converting files:\n","\n".join(files_to_convert)


genome_dict = "None"
converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
converter.create_agp_dict() 

for fname in files_to_convert:
	data = h5dict(fname, mode="r")
	
	#Reading data from dict into numpy arrays
	print "Reading data"
	try:
		strands1 = np.array(data["strands1"],dtype=np.int8)
		cuts1 = np.array(data["cuts1"],dtype=np.uint64)
		chrms1 = np.array(data["chrms1"],dtype=np.uint16)

		strands2 = np.array(data["strands2"],dtype=np.int8)
		cuts2= np.array(data["cuts2"],dtype=np.uint64)
		chrms2 = np.array(data["chrms2"],dtype=np.uint16)

	except:
		raise Exception("Some of the keys were not found in dictionary")
	print "Done with reading"
	
	print "Reading genome"
	
	if genome_dict == "None":
		genome_dict = data["misc"]
		new_genome_dict={"genome":{"label2idx":genome_db_chrmLevel.label2idx,"idx2label":genome_db_chrmLevel.idx2label}}
	elif genome_dict != data["misc"]:
		raise Exception("Possible genome mismatch in file "+fname)
	
	print "Genome mathces!"
	# convert data
	if genome_db_contigLevel.chrmCount >= np.iinfo(np.int16).max - 1:
		raise "Too many chromosomes. Current limit is ",np.iinfo(np.int16).max
		
	for i in xrange(genome_db_contigLevel.chrmCount):
		if i % 100 == 0:
			print i, " done"
		where1 = np.where(chrms1==i)
		where2 = np.where(chrms2==i)
		contigName = "N"+genome_db_contigLevel.idx2label[i]
		if not (contigName in converter.contigsInfo.keys()):
			chrms1[where1] = np.iinfo(np.int16).max
			chrms2[where2] = np.iinfo(np.int16).max
		else:
			chrmName = genome_db_chrmLevel.label2idx[converter.contigsInfo[contigName][0]]
			
			chrms1[where1] = chrmName
			chrms2[where2] = chrmName
			
			if converter.contigsInfo[contigName][3] == "+":
				cuts1[where1] = converter.contigsInfo[contigName][1] + cuts1[where1] - 1 
				cuts2[where2] = converter.contigsInfo[contigName][1] + cuts2[where2] - 1
			elif converter.contigsInfo[contigName][3] == "-":
				cuts1[where1] = converter.contigsInfo[contigName][2] - cuts1[where1] + 1
				cuts2[where2] = converter.contigsInfo[contigName][2] - cuts2[where2] + 1
				strands1[where1] = np.array(np.logical_not(strands1[where1]),dtype=np.int8)
				strands2[where2] = np.array(np.logical_not(strands1[where2]),dtype=np.int8)
			else:
				print converter.contigsInfo[contigName]
				raise
	
	where = np.logical_and(chrms1!=np.iinfo(np.int16).max, chrms2!=np.iinfo(np.int16).max)

	assert sum(where) > len(chrms1)/2 # we lost <= 50% of data
	
	print "Saving"
	out = h5dict(fname+".updGenome", mode="w")
	out.update({"strands1":strands1[where],
				"strands2":strands2[where],
				"cuts1":cuts1[where],
				"cuts2":cuts2[where],
				"chrms1":chrms1[where],
				"chrms2":chrms2[where],
				"misc":new_genome_dict})
	print "Data saved to ",fname