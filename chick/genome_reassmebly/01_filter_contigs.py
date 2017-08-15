from Bio import SeqIO
import os

def checkEnzyme(seq,enzyme="AAGCTT",min_count = 2):
	p = 0
	count = 0
	enzyme_upper = enzyme.upper()
	enzyme_lower = enzyme.lower()
	ind = 0
	while ind<len(seq):
		if seq[ind]==enzyme_upper[p] or seq[ind]==enzyme_lower[p]:
			if p==0:
				ind0 = ind
			p += 1
			if p==len(enzyme):
				p = 0
				count += 1
		else:
			if p != 0:
				ind=ind0
			p = 0
			
		if count == min_count:
			return True
			
		ind += 1
	return False
			
genome_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/"

filtered_contigs = []

totalN = 0
totalLen = 0

addedN = 0
addedLen = 0

for subfolder in ["unplaced_scaffolds/FASTA/","unlocalized_scaffolds/FASTA/"]:
#for subfolder in ["test"]:
	for file in [f for f in os.listdir(genome_folder+"/"+subfolder) if f.split(".")[-1]=="fna"]:
		print "Parsing file ",genome_folder+"/"+subfolder+"/"+file
		records = SeqIO.parse(open(genome_folder+"/"+subfolder+"/"+file),"fasta")
		for ind,record in enumerate(records):
			totalN += 1
			totalLen += len(record.seq)
			if len(record.seq)>50000 and checkEnzyme(record.seq):
				filtered_contigs.append(record)
				addedN+=1
				addedLen+=len(record.seq)

print "-----parsed ",totalN," contigs, ",totalLen," bp"
print "-----passed filters: ",addedN," contigs, ",addedLen," bp"

for subfolder in ["placed_scaffolds/FASTA/"]:
#for subfolder in []:
	for file in [f for f in os.listdir(genome_folder+"/"+subfolder) if f.split(".")[-1]=="fna"]:
		print "Parsing file ",genome_folder+"/"+subfolder+"/"+file
		records = SeqIO.parse(open(genome_folder+"/"+subfolder+"/"+file),"fasta")
		for ind,record in enumerate(records):
			totalN += 1
			totalLen += len(record.seq)
			filtered_contigs.append(record)
			addedN+=1
			addedLen+=len(record.seq)
		
print "-----parsed ",totalN," contigs, ",totalLen," bp"
print "-----passed filters: ",addedN," contigs, ",addedLen," bp"

SeqIO.write(filtered_contigs,genome_folder+"/all_contigs.filtered.fa","fasta")