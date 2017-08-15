#This script based on prepJuiceBoxInput2.py and was designed to use chromosome-based fragment datasets

from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import setExceptionHook
from mirnylib import genome
import numpy as np
import sys
import random
import os
random.seed()
setExceptionHook()

if len(sys.argv) !=4:
	print "Usage:"
	print sys.argv[0]," fragment_dataset_filename enzyme out_file_prefix"
	sys.exit()
else:
	out_file_prefix = sys.argv[3]
	enzyme = sys.argv[2]
	fr_daraset = sys.argv[1]

#Step 1. Prepare chrm.size file

out = open(out_file_prefix+".chrm.size","w")
genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")


for i in xrange(genome_db.chrmCount):
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