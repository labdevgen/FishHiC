import numpy as np
import os
import gzip

import subprocess
from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC
from mirnylib.numutils import fillDiagonal

domain_res='fragment'

genome_name='mm9'
genome_folder='/mnt/storage/home/vsfishman/HiC/fasta/'
genome_db = genome.Genome(genome_folder+genome_name, readChrms=['#', 'X'])

base_filename="Sp_full"
base_folder='/mnt/storage/home/vsfishman/HiC/data/'

maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'

base_out_folder = "/mnt/storage/home/vsfishman/tmp/HiC_tmp/data/"

def Generate_one_chromosome_file(chrNumb):

	o_file=base_out_folder+"fitHiC/i_files/"+base_filename+".fithic"
	fragment_dataset_filename=base_out_folder+"fitHiC/i_files/"+'fragment_dataset_'+base_filename+'_chr'+str(chrNumb)+'.hdf5'

	if not os.path.isfile(fragment_dataset_filename):
		fragments = fragmentHiC.HiCdataset(
			filename=fragment_dataset_filename,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='w')
		fragments.parseInputData(
              dictLike=maped_reads_filepath, removeSS=True)

		fragments.filterRsiteStart(offset=5)
		fragments.filterDuplicates()
		fragments.filterLarge()
		fragments.filterExtreme(cutH=0.005, cutL=0)
	else:
		fragments = fragmentHiC.HiCdataset(
			filename=fragment_dataset_filename,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='a')

	print "Filtering, leaving only chr ", genome_db.idx2label[chrNumb]
	#leave only frgaments from the chrNumb (nterchromosomal)
	fragments.maskFilter((fragments.chrms1 == chrNumb))
	fragments.maskFilter((fragments.chrms2 == chrNumb)) 

	print "Seting RE"
	#Setting info about restriction enzyme, calculating absolute indexes
	fragments.setRfragAbsIdxs('HindIII')
	
	numBins = len(fragments.genome.rsites[chrNumb])
	print "Total numBins (RSites) on chr ",  genome_db.idx2label[chrNumb], " = ", numBins
	
	rfragAbsIdxs1 = fragments.rfragAbsIdxs1-fragments.genome.chrmStartsRfragCont[chrNumb]
	rfragAbsIdxs2 = fragments.rfragAbsIdxs2-fragments.genome.chrmStartsRfragCont[chrNumb]
	print "Total number of fragments = ",len(rfragAbsIdxs1)
	
	if len(rfragAbsIdxs1) != len(rfragAbsIdxs2):
		print "rfragAbsIdxs1=",rfragAbsIdxs1
		print "rfragAbsIdxs2=",rfragAbsIdxs2
		print "len(rfragAbsIdxs1)=", len(rfragAbsIdxs1)
		print "len(rfragAbsIdxs2)=", len(rfragAbsIdxs2)	
		raise "FRAGMENT INDEXING ERROR 1!!!"
	if (min(rfragAbsIdxs1)<0 or min(rfragAbsIdxs2)<0 ):
		print "min(rfragAbsIdxs1)=",min(rfragAbsIdxs1)
		print "min(rfragAbsIdxs2)=",min(rfragAbsIdxs2)	
		raise "FRAGMENT INDEXING ERROR 2!!!"
	if (max(rfragAbsIdxs1)> numBins-1 or max(rfragAbsIdxs2) > numBins-1):
		print "max (rfragAbsIdxs1)=", max (rfragAbsIdxs1)
		print "max (rfragAbsIdxs2)=", max (rfragAbsIdxs2)
		print "numBins=",numBins
		raise "FRAGMENT INDEXING ERROR 3!!!"

	print "FRAGMENT INDEXING - passed"

	#Creating label array
	label = np.array(rfragAbsIdxs1,dtype='int64')
	label *= numBins
	label += rfragAbsIdxs2

	#Creating count array
	counts = np.bincount(label, minlength=numBins ** 2)
	counts.shape = (numBins, numBins)

	#Counting
	for i in xrange(len(counts)):
		counts[i, i:] += counts[i:, i]
		counts[i:, i] = counts[i, i:]

	#Filling diagonal reads
	
	#diag = np.diag(counts)
	#fillDiagonal(counts, diag/2)
	fillDiagonal(counts, 0)

	BinsToDescribe=np.zeros(numBins) # Info about which RSites should be described in .fragments file later

#	f_out = gzip.open (o_file+"_chr"+str(chrNumb)+".contacts.zip","w")	
	f_out = open (o_file+"_chr"+str(chrNumb)+".contacts.zip","w")	
	print "Writing file ",o_file+"_chr"+str(chrNumb)+".contacts.zip"
	for i in range(numBins-1):
		for j in range(i+1,numBins):
			if (counts[i,j] != 0):
				s = ""
				s += str(chrNumb) + "\t"
				s += str(fragments.genome.rfragMids[chrNumb][i]) + "\t"
				s += str(chrNumb) + "\t"
				s += str(fragments.genome.rfragMids[chrNumb][j]) + "\t"
				s += str(counts[i,j])+"\n"
				f_out.write(s)
				BinsToDescribe[i] = 1
				BinsToDescribe[j] = 1				

	f_out.close()
	
#	f_out = gzip.open (o_file+"_chr"+str(chrNumb)+".fragments.zip","w")
	f_out = open (o_file+"_chr"+str(chrNumb)+".fragments.zip","w")
	print "Writing file ",o_file+"_chr"+str(chrNumb)+".fragments.zip"
	
	for ind,val in enumerate(BinsToDescribe):
		if (val == 1):
			s = ""
			s += str(chrNumb) + "\t0\t"		
			s += str(fragments.genome.rfragMids[chrNumb][ind])+"\t"
			s += str(sum(counts[ind]))+"\t"
			s += "1\n"
			f_out.write(s)
	f_out.close()

def Generate_whole_genome_chromosome_file(mappability = 1): #TODO - use mappability, now it's always =1
	fragment_dataset_filename=base_folder+'fragment_dataset_'+base_filename+'.hdf5'
	o_file = base_out_folder + "fitHiC/i_files/"+base_filename

	if not os.path.isfile(fragment_dataset_filename):
		fragments = fragmentHiC.HiCdataset(
			filename=fragment_dataset_filename,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='w')
		fragments.parseInputData(
              dictLike=maped_reads_filepath, removeSS=True)

		fragments.filterRsiteStart(offset=5)
		fragments.filterDuplicates()
		fragments.filterLarge()
		fragments.filterExtreme(cutH=0.005, cutL=0)
	else:
		print "fragment_dataset "+ fragment_dataset_filename + " found.\n IMPORTANT: considering all requerd filters have been done"
		fragments = fragmentHiC.HiCdataset(
			filename=fragment_dataset_filename,
			genome=genome_db,
			maximumMoleculeLength=500,
			mode='a')

	print "Seting RE"
	#Setting info about restriction enzyme, calculating absolute indexes
	fragments.setRfragAbsIdxs('HindIII')

	rfragAbsIdxs1 = fragments.rfragAbsIdxs1
	rfragAbsIdxs2 = fragments.rfragAbsIdxs2
	print "Total number of fragments = ",len(rfragAbsIdxs1)

	if len(rfragAbsIdxs1) != len(rfragAbsIdxs2):
		print "rfragAbsIdxs1=",rfragAbsIdxs1
		print "rfragAbsIdxs2=",rfragAbsIdxs2
		print "len(rfragAbsIdxs1)=", len(rfragAbsIdxs1)
		print "len(rfragAbsIdxs2)=", len(rfragAbsIdxs2)	
		raise "FRAGMENT INDEXING ERROR 1!!!"
	if (min(rfragAbsIdxs1)<0 or min(rfragAbsIdxs2)<0 ):
		print "min(rfragAbsIdxs1)=",min(rfragAbsIdxs1)
		print "min(rfragAbsIdxs2)=",min(rfragAbsIdxs2)	
		raise "FRAGMENT INDEXING ERROR 2!!!"
	if (max(rfragAbsIdxs1)> fragments.genome.chrmEndsRfragCont[-1] or max(rfragAbsIdxs2) > fragments.genome.chrmEndsRfragCont[-1]):
		print "max (rfragAbsIdxs1)=", max (rfragAbsIdxs1)
		print "max (rfragAbsIdxs2)=", max (rfragAbsIdxs2)
		print "numBins=",fragments.genome.chrmEndsRfragCont[-1]
		raise "FRAGMENT INDEXING ERROR 3!!!"

	print "FRAGMENT INDEXING - passed"

	print "Initialyzing heatmap"
	max_rsites_number = max([len(i) for i in fragments.genome.rsites])
	heatmap=np.zeros(shape=(fragments.genome.chrmCount,max_rsites_number,max_rsites_number),dtype=np.uint16)
	print "Max numBins (RSites) in chromosome= ", max_rsites_number
	


	print "Creating chrRsites array"
	RsiteToChr=np.zeros(max_rsites_number*fragments.genome.chrmCount,dtype=np.uint8)
	for i in range(0,fragments.genome.chrmCount):
		RsiteToChr[fragments.genome.chrmStartsRfragCont[i]:fragments.genome.chrmEndsRfragCont[i]] = i
	
	print "Filling heatmap"
	BinsToDescribe=np.zeros(shape=(fragments.genome.chrmCount,max_rsites_number),dtype=np.int8)	

	l = len(rfragAbsIdxs1)
	for i in xrange(l):
		if (i % (l/10)) == 0:
			print (i / (l/10)),"0 %"
			
		ChrN1=RsiteToChr[rfragAbsIdxs1[i]]
		ChrN2=RsiteToChr[rfragAbsIdxs2[i]]
		if ( ChrN1 == ChrN2): #if it is intrachromosomal contact
			rfragAbsIdxs1_onChr=rfragAbsIdxs1[i]-fragments.genome.chrmStartsRfragCont[ChrN1]
			rfragAbsIdxs2_onChr=rfragAbsIdxs2[i]-fragments.genome.chrmStartsRfragCont[ChrN1]			
			BinsToDescribe[ChrN1][rfragAbsIdxs1_onChr] = 1
			BinsToDescribe[ChrN2][rfragAbsIdxs2_onChr] = 1
			heatmap[ChrN1][rfragAbsIdxs1_onChr][rfragAbsIdxs2_onChr] += 1
			if heatmap[ChrN1][rfragAbsIdxs1_onChr][rfragAbsIdxs2_onChr] >= 64000:
				raise "Type int16 used in heatmap is not compatible with N of contact >64000"
	
	
	f_out = open (o_file+".contacts","w")
	
	
	print "Total number of non-empty bins (rfrags)=",np.sum(BinsToDescribe)

	print "Writing file ",o_file+".contacts"
	f_out = open (o_file+".contacts","w")
	for i in xrange(l):
		if (i % (l/10)) == 0:
			print (i / (l/10)),"0 %"
			
		ChrN1=RsiteToChr[rfragAbsIdxs1[i]]
		ChrN2=RsiteToChr[rfragAbsIdxs2[i]]

		if ( ChrN1 == ChrN2): #if it is intrachromosomal contact
			rfragAbsIdxs1_onChr=rfragAbsIdxs1[i]-fragments.genome.chrmStartsRfragCont[ChrN1]
			rfragAbsIdxs2_onChr=rfragAbsIdxs2[i]-fragments.genome.chrmStartsRfragCont[ChrN1]			
			if (heatmap[ChrN1][rfragAbsIdxs1_onChr][rfragAbsIdxs2_onChr] != -1):
				s = ""
				s += str(i) + "\t"
				s += str(fragments.genome.rfragMids[ChrN1][rfragAbsIdxs1_onChr]) + "\t"
				s += str(i) + "\t"
				s += str(fragments.genome.rfragMids[ChrN2][rfragAbsIdxs2_onChr]) + "\t"
				s += str(heatmap[ChrN1,rfragAbsIdxs1_onChr,rfragAbsIdxs2_onChr])+"\n"
				heatmap[ChrN1][rfragAbsIdxs1_onChr][rfragAbsIdxs2_onChr] = -1				
				f_out.write(s)

				
	f_out.close()
	
	f_out = open (o_file+".fragments","w")

	print "Writing file ",o_file+".fragments"
	for i in range(fragments.genome.chrmCount):	
		for j in xrange(max_rsites_number):	
			if (BinsToDescribe[i][j] == 1):
				s = ""
				chrNumb = i
				s += str(chrNumb) + "\t0\t"		
				s += str(fragments.genome.rfragMids[i][j])+"\t"
				s += str(sum(heatmap[i][j]))+"\t"
				s += "1\n"
				f_out.write(s)
	f_out.close()


def executeBashCommand(bashCommand,doPrint=True):
  print "executing command %s \n" % (bashCommand)
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  if doPrint:
    print p.communicate()[0]
    return ""
  else:
    return p.communicate()[0]
	
def doFitHiC(chrNumb):
	o_file=base_folder+"fitHiC/i_files/"+base_filename+".fithic"
	o_file+= "_chr"+str(chrNumb)
	c = "python2.7 /mnt/storage/home/vsfishman/HiC/fit-hi-c/bin/fit-hi-c.py -f "+o_file+".fragments.zip -i "+o_file+".contacts.zip -o "+base_folder+"fitHiC/o_files/"+" -l"+base_filename+"_"+str(chrNumb)+" -U 11000000"
	print "Executing command ",c
	executeBashCommand(c)

def doFitHiC_whole_genome():
	o_file = base_out_folder + "fitHiC/i_files/"+base_filename
	c = "python2.7 /mnt/storage/home/vsfishman/HiC/fit-hi-c/bin/fit-hi-c.py -f "+o_file+".fragments.zip -i "+o_file+".contacts.zip -o "+base_out_folder+"fitHiC/o_files/"+" -l"+base_filename+" -U 5000000"
	print "Executing command ",c
	executeBashCommand(c)

	
#base_names=["Sp_full","ESC_full","Fib_full"]
base_names=["Sp_1000","ESC_full","Fib_full"]
chrNumbs=[0]
for i in chrNumbs:
	for j in base_names:
		base_filename=j
		Generate_whole_genome_chromosome_file()
		doFitHiC_whole_genome()
#		Generate_one_chromosome_file(i)
#		doFitHiC(i)