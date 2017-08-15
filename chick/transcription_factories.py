import matplotlib
matplotlib.use('Agg')
import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome
from hiclib import fragmentHiC 
from mirnylib import plotting
from hiclib import binnedData
logging.basicConfig(level=logging.DEBUG)

import matplotlib.pyplot as plt
import numpy as np

from mirnylib.systemutils import setExceptionHook
setExceptionHook()
import datetime

#######################


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

data = {
	"FibsTop" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII_refined.frag",
				"/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.Fibs.testFPKM.ChrmLevel.regions"],
#	"BloodTop" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII_refined.frag",
#			"/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.Blood.topFPKM.ChrmLevel.regions"],
	"FibsLow" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII_refined.frag",
				"/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.Fibs.test2FPKM.ChrmLevel.regions"],
#	"BloodLow" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII_refined.frag",
#			"/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.Blood.lowFPKM.ChrmLevel.regions"]

}


###########################

for tmp_ind,key in enumerate(sorted(data)):
	if (tmp_ind == 0) or data[sorted(data)[tmp_ind-1]] != key:
		fragments = fragmentHiC.HiCdataset(
			filename=data[key][0]+".tmp",
			inMemory=True,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='w',
			enzymeName="HindIII")
		fragments.load(data[key][0])

		assert fragments._isSorted()
		fragids1 = fragments._getVector("fragids1")
		fragids2 = fragments._getVector("fragids2")
		fs  = fragments.fragmentSum() #contains total N of contacts for each fragment

	regions = np.genfromtxt(data[key][1],dtype=None)
	desierd_fragids = [] #contains ids of fragment froma desierd regions
						# id is a nt of rfrag midpoint + chrm*genome_db.fragIDmult
	print "---------------------------Processing ",key,"--------------------"
	print "defining desired fragments"
	for region in regions:
		chrm,st,end = region
		if not chrm in genome_db.label2idx.keys():
			print chrm
			continue
		else:
			chrm = genome_db.label2idx[chrm]
		
		if st == genome_db.chrmLens[chrm]:
			rFrgagChromIdx = len(genome_db.rsites[chrm])-1
		else:
			rFrgagChromIdx = np.searchsorted(genome_db.rsites[chrm],st)
		rFrgagContIdx_left = genome_db.rfragMids[chrm][rFrgagChromIdx] + chrm*genome_db.fragIDmult
		assert rFrgagContIdx_left in fragments.rFragIDs
		
		if end >= genome_db.chrmLens[chrm]:
			rFrgagChromIdx = len(genome_db.rsites[chrm])-1
		else:
			rFrgagChromIdx = np.searchsorted(genome_db.rsites[chrm],end)
		rFrgagContIdx_riht = genome_db.rfragMids[chrm][rFrgagChromIdx] + chrm*genome_db.fragIDmult
		assert rFrgagContIdx_riht in fragments.rFragIDs

		if rFrgagContIdx_left==rFrgagContIdx_riht:
			desierd_fragids.append(rFrgagContIdx_riht)
		else:
			fr = np.where(fragments.rFragIDs==rFrgagContIdx_left)[0][0]
			to = np.where(fragments.rFragIDs==rFrgagContIdx_riht)[0][0]
#			desierd_fragids += list(fragments.rFragIDs[fr:to])
			desierd_fragids.append(fragments.rFragIDs[fr])
			
	desierd_fragids = np.unique(desierd_fragids)
	print "Defined ",len(desierd_fragids)," desierd fragids"
	print desierd_fragids

	print "Fragments defined. Filtering"

	#desierd_frags - bool array to index fragids1 and fragids2
	#array contain Trus if fragids1 and fragids2 both are in desierd_fragids are intrachromosomal
	desierd_frags = np.logical_and(np.in1d(fragids1,desierd_fragids),np.in1d(fragids2,desierd_fragids)) 
#	desierd_frags = np.logical_and(desierd_frags,np.not_equal(fragments._getVector("chrms1"),fragments._getVector("chrms2")))
	desired_indexes = np.nonzero(desierd_frags) #indexes of elements(contacts) in fragids1 and fragids2, that we are goint to proceed
	print "Filtering done"
	print len(desired_indexes[0])," contacts to process"	
	
	print "Calculating weights"
	
	assert len(fs) == len(fragments.rFragIDs)


	
	#weight of fragment will be calulated as 1/(s1*s2) where s1 and s2 are total sum of contacts for this fragment
	#fragments._getVector("rfragAbsIdxs1") will return array with the same lengths as fragids1 or 2
	#each element of rfragAbsIdxs1 contains not an id of fragment (recalling, id is midpoint + chrm*genome_db.fragIDmult),
	#but absolute index, so that 1st fragment has an index 0 (whereas id is ~2600), second has index 1(id ~4500) and etc.
	#then subset rfragAbsIdxs1/2 using desired_indexes to get only indexes of desierd_frags
	#finally, get total N of contacts fot there fragments indexing fs array
	w1 = fs[fragments._getVector("rfragAbsIdxs1")[desired_indexes]]
	w2 = fs[fragments._getVector("rfragAbsIdxs2")[desired_indexes]]
	w = 1. / np.multiply(w1,w2)

	print "Calculating normalizaed contacts"
	print key,np.average(w),np.std(w),min(w),max(w),np.median(w),sum(w)/len(desired_indexes)