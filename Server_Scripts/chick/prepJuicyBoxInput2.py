from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import setExceptionHook
from mirnylib import genome
import numpy as np
import sys
import random
import os
random.seed()
setExceptionHook()

if not (len(sys.argv) in [4,5]):
	print "Usage:"
	print sys.argv[0]," fragment_dataset_filename enzyme out_file_prefix mode(optional)"
	print sys.argv[0],"mode can be eaither 'chr' or empty"
	sys.exit()
else:
	out_file_prefix = sys.argv[3]
	enzyme = sys.argv[2]
	fr_daraset = sys.argv[1]
if len(sys.argv) == 5:
	mode = sys.argv[4]
	if mode != "chr":
		raise Exception("Unknown mode "+mode+"\nMode can only be chr")
	agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP"
else:
	mode = "all"

#Step 1. Load agp info about genome assmebly, for "chr" mode only
#For "all" mode Load genome 
#and prepare chrm.size file

out = open(out_file_prefix+".chrm.size","w")
genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")

# genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				# readChrms=[],
				# chrmFileTemplate="%s.fna")


if mode == "chr":
	agp_dict = {}
	agp_files = [f for f in os.listdir(agp_folder) if f.endswith(".agp")]
	for f in sorted(agp_files,key=len):
		with open(agp_folder+"/"+f) as agp_file:
			chr_len = 0
			for line in agp_file:
				if line[0]=="#": continue
				line = line.split()
				if line[4] == "O" or line[4] == "W": #O is a contig, U is a gap, N is a centromere, W - contig on chrm W
					if line[-1] == "+":
						agp_dict[line[5]] = [f.split(".")[0],int(line[1]),line[-1]] #contigName --> [chr,contigStartBp,orientataion]
					elif line[-1] == "-":
						agp_dict[line[5]] = [f.split(".")[0],int(line[2]),line[-1]] #contigName --> [chr,contigEndBp,orientataion]
					else:
						raise
				elif line[4] != "U" and line[4] != "N":
					raise
				assert chr_len < int(line[2])
				chr_len = int(line[2])
			out.write(f.split(".")[0]+"\t"+str(chr_len)+"\n")
else:	
	for i in xrange(genome_db.chrmCount):
		if genome_db.chrmCount > 100:
			out.write("N"+genome_db.idx2label[i]+"\t"+str(genome_db.chrmLens[i])+"\n")
		else:
			out.write(genome_db.idx2label[i]+"\t"+str(genome_db.chrmLens[i])+"\n")
out.close()

#Step 2. Load fragments and save it in juicybox format

Fr = HiCdataset(fr_daraset+".tmp"+str(random.randint(0,10000)),enzymeName = enzyme,
                       mode='w', genome=genome_db,inMemory=True)
Fr.load(fr_daraset)

out=open(out_file_prefix+".contacts","w")
print "Reporting ",Fr.N," contacts"
dotsize = Fr.N / 500
print "one dot is ",dotsize," contacts"

strands1 = np.array(Fr._getVector("strands1"),dtype=np.int8)
cuts1 = np.array(Fr._getVector("cuts1"),dtype=np.uint64)
chrms1 = np.array(Fr._getVector("chrms1"),dtype=np.uint16)
fragids1 = np.array(Fr._getVector("rfragAbsIdxs1"),dtype=np.uint32)

strands2 = np.array(Fr._getVector("strands2"),dtype=np.int8)
cuts2= np.array(Fr._getVector("cuts2"),dtype=np.uint64)
chrms2 = np.array(Fr._getVector("chrms2"),dtype=np.uint16)
fragids2 = np.array(Fr._getVector("rfragAbsIdxs2"),dtype=np.uint32)

#for debug use only----------------
#check_idxs = [random.randint(0,len(strands1)-1) for i in xrange(20)]
#strands1c,strands2c,cuts1c,cuts2c,chrms1c,chrms2c = strands1[check_idxs],strands2[check_idxs],cuts1[check_idxs],cuts2[check_idxs],chrms1[check_idxs],chrms2[check_idxs]
#----------------

if mode == "chr":
	print "Converting contigs to assembled chromosomes"
	#prepare list of not assembled contigs
	unassembled_list = []
	chrmLables = sorted(agp_dict.keys())
	for i in xrange(genome_db.chrmCount):
		if not ("N"+genome_db.idx2label[i] in chrmLables):
			unassembled_list.append(i)
	
	assert len(unassembled_list)<genome_db.chrmCount
	
	sample = chrms1[::len(chrms1)/5000] #random sample of chrm1
	Nremoved = sum([i in unassembled_list for i in sample])
	Premoved = 1.-(1.-float(Nremoved)/len(sample))**2
	print "Going to remove ~",Premoved*100,"% of all contacts" #1-q^2 - probability of at least one success in 2 tests (one for chrm1 and one for chrm2) 
	assert (Premoved < 0.5)
	

	print "Creating new chrms arrays"
	chrms_new_1 = np.array(chrms1,dtype="|S8")
	chrms_new_2 = np.array(chrms2,dtype="|S8")	
	print "Done"
	
	count = 0
	for i in xrange(genome_db.chrmCount):
		if i % 100 == 0:
			print i, " done"
		where1 = np.where(chrms1==i)
		where2 = np.where(chrms2==i)
		contigName = "N"+genome_db.idx2label[i]
		if not (contigName in chrmLables):
			chrms_new_1[where1] = "UN"
			chrms_new_2[where2] = "UN"
		else:
			chrmName = agp_dict[contigName][0]
			
			chrms_new_1[where1] = chrmName
			chrms_new_2[where2] = chrmName
			
			if agp_dict[contigName][2] == "+":
				cuts1[where1] = agp_dict[contigName][1] + cuts1[where1] - 1 
				cuts2[where2] = agp_dict[contigName][1] + cuts2[where2] - 1				
			else:
				cuts1[where1] = agp_dict[contigName][1] - cuts1[where1] + 1
				cuts2[where2] = agp_dict[contigName][1] - cuts2[where2] + 1
				strands1[where1] = np.array(np.logical_not(strands1[where1]),dtype=np.int8)
				strands2[where2] = np.array(np.logical_not(strands1[where2]),dtype=np.int8)
				
	chrms1 = chrms_new_1
	chrms2 = chrms_new_2
	
	#DEBUG
#	print "-----DEBUG----"
#	for indc,ind in enumerate(check_idxs):
#		print "\t".join(map(str,[strands1c[indc],strands2c[indc],cuts1c[indc],cuts2c[indc],genome_db.idx2label[chrms1c[indc]],genome_db.idx2label[chrms2c[indc]]]))
#		print "\t".join(map(str,[strands1[ind],strands2[ind],cuts1[ind],cuts2[ind],chrms1[ind],chrms2[ind]]))
#		print "-----"
	#END DEBUG	
print "Sorting array"

#sorted_index = np.lexsort((chrms2,chrms1))
#np.savetxt("chrms1.txt",chrms1)
#np.savetxt("chrms2.txt",chrms2)
sorted_index = np.lexsort((np.where(chrms1>chrms2,chrms1,chrms2),np.where(chrms1<=chrms2,chrms1,chrms2)))
print "Done"

print "Start reporting"

count = 0

for i in sorted_index:
	count += 1			
	if count % dotsize == 0:
		sys.stdout.write(".")
		sys.stdout.flush()
#		if count % (dotsize*10) == 0:
#			print "\n",count," processed"
	if mode == "chr":
		if chrms1[i] == "UN":
			continue
		if chrms2[i] == "UN":
			continue
		s = str(i) + "\t" +\
			str(strands1[i])+"\t"+\
			chrms1[i] + "\t" + \
			str(cuts1[i])+"\t"+ \
			str(fragids1[i])+"\t"+\
			str(strands2[i])+"\t"+\
			chrms2[i] + "\t" +\
			str(cuts2[i]) + "\t   " +\
			str(fragids2[i]) + "\t" +\
			"0\t0\n"
	else:
		s = str(i) + "\t" +\
			str(strands1[i])+"\t"+\
			"N"+genome_db.idx2label[chrms1[i]] + "\t" + \
			str(cuts1[i])+"\t"+ \
			str(fragids1[i])+"\t"+\
			str(strands2[i])+"\t"+\
			"N"+genome_db.idx2label[chrms2[i]] + "\t" +\
			str(cuts2[i]) + "\t" +\
			str(fragids2[i]) + "\t" +\
			"0\t0\n"
	out.write(s)

out.close()