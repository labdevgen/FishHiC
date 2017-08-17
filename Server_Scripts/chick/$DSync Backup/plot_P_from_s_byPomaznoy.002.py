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
import scipy
linr=scipy.stats.linregress

from mirnylib.systemutils import setExceptionHook
setExceptionHook()
import datetime

#######################


genome_db_contig = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/",
				readChrms=[],
				chrmFileTemplate="N%s.fa")

genome_db_chrm = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")


data = {
#"FibsChrmLevel" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII_refined.frag", genome_db_chrm],
"Fibs_micro" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII_refined.frag", 
					genome_db_contig,"micro"],
"Fibs_macro" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII_refined.frag", 
					genome_db_contig,"macro"],
"Fibs_all" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII_refined.frag", 
					genome_db_contig,None],
#"Fibs_active" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/ChEF-all-HindIII_refined.frag", 
#					genome_db_contig,"/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.Fibs.topFPKM.regions"],
"Blood_all" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII_refined.frag", 
					genome_db_contig,None],
#"Blood_active" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII_refined.frag", 
#					genome_db_contig,"/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff.Blood.topFPKM.regions"],
"Blood_micro" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII_refined.frag", 
					genome_db_contig,"micro"],
"Blood_macro" : ["/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filtered/Blood-all-HindIII_refined.frag", 
					genome_db_contig,"macro"]
}


###########################

abscissa={}
ordinata={}

for key in data:
	genome_db = data[key][1]
	fragments = fragmentHiC.HiCdataset(
			filename=data[key][0]+".tmp",
			inMemory=True,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='w',
			enzymeName="HindIII")
	fragments.load(data[key][0])
	print "Preparing regions for ",key
	
	chrms1 = np.array(fragments._getVector("chrms1"),dtype=np.uint16)
	chrms1 = chrms1[::len(chrms1)/10000] #This is a representative sample of 0.01% of contacts
	if data[key][2] == None:  #data[key][2] possible values: micro, macro, None
								#None counts all chromosomes
		regions = []
		for chrm in range(genome_db.chrmCount):
			if sum(chrms1==chrm) > 1: #Check that there are enough contacts on this chromsome, 
										#i.e. contacts of this chromosome account for > 0.01% of all contacts
				regions.append((chrm,0,genome_db.chrmLens[chrm]))
		(abscissa[key],ordinata[key])=fragments.plotScaling(withinArms=False,regions=regions)
	elif (data[key][2] == "micro") or (data[key][2] == "macro"):
		regions = []
		from lib_coordinates_converter import *
		converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
		converter.prepare_acc2chrmDict()
		converter.create_agp_dict()
		chrmLabels = converter.getChrmsByType(format="contig")[data[key][2]]
		print "Using chromosomes: ",chrmLabels
		chrms = [genome_db.label2idx[c[1:]] for c in chrmLabels if c[1:] in genome_db.label2idx.keys()] #Stupid names: in genome_db they are T_nnnn due to chrmFileTemplate="N%s.fa", but normally they are NT_nnnn.
		print len(chrms)," chrms out of ",len(chrmLabels)," found in genome"
		for chrm in range(genome_db.chrmCount):
			if sum(chrms1==chrm) > 1 and (chrm in chrms):
				regions.append((chrm,0,genome_db.chrmLens[chrm]))
		print "After filtering using ",len(regions)," regions"
		(abscissa[key],ordinata[key])=fragments.plotScaling(withinArms=False,regions=regions)
	elif (data[key][2]=="chr17_chr19_chr_22"):
		regions = []
		from lib_coordinates_converter import *
		converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
		converter.prepare_acc2chrmDict()
		converter.create_agp_dict()
		chrmLabels = [i for chr in ["chr17","chr19","chr22"] for i in converter.chrm2contig[chr]] 
		print "Using chromosomes: ",chrmLabels
		chrms = [genome_db.label2idx[c[1:]] for c in chrmLabels if c[1:] in genome_db.label2idx.keys()] #Stupid names: in genome_db they are T_nnnn due to chrmFileTemplate="N%s.fa", but normally they are NT_nnnn.
		print len(chrms)," chrms out of ",len(chrmLabels)," found in genome"
		for chrm in range(genome_db.chrmCount):
			if sum(chrms1==chrm) > 1 and (chrm in chrms):
				regions.append((chrm,0,genome_db.chrmLens[chrm]))
		print "After filtering using ",len(regions)," regions"
		(abscissa[key],ordinata[key])=fragments.plotScaling(withinArms=False,regions=regions)

	else: 
		from hiclib.hicShared import h5dictBinarySearch
		assert fragments._isSorted()
		regions = np.genfromtxt(data[key][2],dtype=None)
		all_chrms = np.unique(regions["f0"])
		bad_chrms = []
		desierd_fragids = []
		
		c1_h5 = fragments.h5dict.get_dataset("chrms1")
		p1_h5 = fragments.h5dict.get_dataset("cuts1")
		c2_h5 = fragments.h5dict.get_dataset("chrms2")
		p2_h5 = fragments.h5dict.get_dataset("cuts2")

		for ind,region in enumerate(regions):
			chrm,st,end = region
			chrm = chrm[1:] #Stupid names: in genome_db they are T_nnnn due to chrmFileTemplate="N%s.fa", but normally they are NT_nnnn.
			if not chrm in genome_db.label2idx.keys():
				continue
			else:
				chrm = genome_db.label2idx[chrm]
				low1 = h5dictBinarySearch(c1_h5, p1_h5, (chrm, st), "left")
				high1 = h5dictBinarySearch(c1_h5, p1_h5, (chrm, end), "right")
				if low1 != high1:
					desierd_fragids += list(np.unique(fragments._getVector("fragids1", low1, high1)))

				low2 = h5dictBinarySearch(c2_h5, p2_h5, (chrm, st), "left")
				high2 = h5dictBinarySearch(c2_h5, p2_h5, (chrm, end), "right")
				if low2 != high2:
					desierd_fragids += list(np.unique(fragments._getVector("fragids1", low2, high2)))
				
		desierd_fragids = np.unique(desierd_fragids)
		print len(desierd_fragids) #DEBUG
		(abscissa[key],ordinata[key])=fragments.plotScaling(withinArms=False,fragids1=None,fragids2=desierd_fragids)
plt.clf()

out_file = open("P_s.log.txt","a")
out_file.write("Start\t"+str(datetime.datetime.now())+"\n")
for key in data:
	out_file.write(key+"\n")
	out_file.write(str(data[key])+"\n")
	out_file.write(str(abscissa[key])+"\n")
	out_file.write(str(ordinata[key])+"\n")

for key in data:
	ind1=abscissa[key]<7000000
	ind2=abscissa[key]>500000
	ind=ind1*ind2
	slope, intercept, r_value, p_value, std_err = linr(np.log10(abscissa[key][ind]),np.log10(ordinata[key][ind]))
	plt.plot(abscissa[key],ordinata[key], label=key+'	~s**'+str(slope), linewidth=2)

ax = plt.gca()
plt.xscale('log')
plt.yscale('log')

fs = 12
plt.xlabel("Genomic distance (bases)", fontsize=12)
plt.ylabel("Contact probability", fontsize=12)
for xlabel_i in ax.get_xticklabels():
				xlabel_i.set_fontsize(fs)
for xlabel_i in ax.get_yticklabels():
				xlabel_i.set_fontsize(fs)
				
				
xs=abscissa[abscissa.keys()[0]]
ys=1/xs*0.5
plt.plot(xs,ys,label="1/s",linewidth=4)

				
legend = plt.legend(loc=0, prop={"size": 12})
legend.draw_frame(False)
 
plt.savefig('P_from_s.png',dpi=300)
		
