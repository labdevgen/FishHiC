import numpy as np
import sys
import random
random.seed()
from Bio import SeqIO
from lib_coordinates_converter import *


from mirnylib.systemutils import setExceptionHook
setExceptionHook()


def generate_domains_from_file(domains_file):
	domains = np.genfromtxt(domains_file,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
	domains = np.sort(domains,order=["chrm","start"])
	return domains


def get_fasta_fname(fname,fastas_path= "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/galGal5_all_contigs.filtered/"
):
	for i in os.listdir(fastas_path):
		if fname.upper() in i.upper():
			return fastas_path+i
	return -1

converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")

converter.prepare_acc2chrmDict()							
converter.create_agp_dict()						

domains_file = sys.argv[1]
domains = list(generate_domains_from_file(domains_file))
borders = []
for i in xrange(1,len(domains)):
	if domains[i][0]==domains[i-1][0]:
		borders.append((domains[i][0],domains[i-1][2],domains[i][1]))
		assert borders[-1][1]<=borders[-1][2]

N_domains = len(domains)
print len(domains)," defined"

ctotal = 0
cDifChr = 0
total_size = 0
f_out = open(domains_file+".borders_seq.fa","w")
while len(borders) > 0 and total_size <= 40000000:
	d = random.randint(0,len(borders)-1)
	chr,st,end = borders[d]
	del borders[d]
	ctotal += 1
	st = converter.chrmTOcontig(chr,st+1)
	end = converter.chrmTOcontig(chr,end-1)
	
	if (st==None) or (end==None):
		cDifChr += 1
		continue
		
	if st.chrm != end.chrm:
		cDifChr += 1
		continue
	if st.nt > end.nt:
		temp = st
		st = end
		end = temp
	fname = get_fasta_fname(st.chrm)
	assert fname != -1
	records = list(SeqIO.parse(fname,'fasta'))
	assert len(records)==1
	new_st = max(0,st.nt-20000)
	new_end = min(len(records[0].seq),end.nt+20000)
	f_out.write(">"+str(ctotal)+"\n")
	f_out.write(str(records[0].seq[new_st:new_end])+"\n")
	total_size += new_end-new_st

f_out.close()		
print ctotal," domains out of ",N_domains," processed, ",cDifChr," excluded due to border located on different contigs"
print "lenth of border: ",total_size
