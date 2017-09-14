import datetime

now = datetime.datetime.now()

print "Starting",str(now)
import os
from mirnylib.systemutils import setExceptionHook
setExceptionHook()

from mirnylib import genome
from hiclib import binnedData
from mirnylib import h5dict
import numpy as np
import argparse

#genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
#				readChrms=[],
#				chrmFileTemplate="N%s.fa")

genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				 readChrms=[],
				 chrmFileTemplate="%s.fna")


parser = argparse.ArgumentParser()
parser.add_argument("--hmap")
parser.add_argument("--noexp",action='store_true')
args = parser.parse_args()


#domain_res=args.res
hmap=str(args.hmap).strip()

#trying to get resolution from hmap
raw_heatmap = h5dict.h5dict(hmap, mode='r') 
res = int(raw_heatmap['resolution'])
print "resolution defined by heatmap: ",res

BD = binnedData.binnedData(res, genome_db)
print datetime.datetime.now()," loading hmap"
BD.simpleLoad(hmap, 'heatmap')

print datetime.datetime.now()," filtering"
BD.removeBySequencedCount(0.5)
# Remove 1% of regions with low coverage.
BD.removePoorRegions(cutoff=1)
# Truncate top 0.05% of interchromosomal counts (possibly, PCR blowouts).
BD.truncTrans(high=0.0005)

#Remove diagonal
BD.removeDiagonal()

print datetime.datetime.now()," performing IC"
# Perform iterative correction.
BD.iterativeCorrectWithoutSS()

#Export hmap
print datetime.datetime.now()," exporting"
BD.export("heatmap",hmap+".IC")

if not args.noexp:
	if not os.path.exists(hmap+".gzipped_matrix"):
		os.makedirs(hmap+".gzipped_matrix")

	print datetime.datetime.now()," saving as txt"	
	for chrm,start,end in zip(range(len(genome_db.chrmStartsBinCont)),genome_db.chrmStartsBinCont,genome_db.chrmEndsBinCont):
		np.savetxt(hmap+".gzipped_matrix/N"+genome_db.idx2label[chrm]+".gz",BD.dataDict["heatmap"][start:end,start:end],delimiter="\t")
print datetime.datetime.now()," Done."