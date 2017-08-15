import os
from mirnylib.genome import Genome 


genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")


out = open(genome_db.genomePath+"/GalGal5ChrmLevel.fai","w")

for ind,faiFail in enumerate(genome_db.fastaNames):
	try:
		lines = open(faiFail+".fai").readlines()
	except:
		raise Exception("Fai file not found. Run samtools faidx on ",faiFail," first")
	if len(lines)>1:
		print faiFail,lines
		raise
	l = lines[0].strip().split()
	l[0] = str(ind)
	out.write("\t".join(l)+"\n")
	
out.close()