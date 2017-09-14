import sys
import numpy as np
from intervaltree import Interval, IntervalTree

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
import figPath
figure_path=figPath.figure_path+"FISH_domains.txt"

from mirnylib.systemutils import setExceptionHook
setExceptionHook()

from lib_coordinates_converter import *
converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
converter.create_agp_dict() 


class TAD():
	def __init__(self,chr,start,end,type="NA"):
		self.chr = chr
		self.start = start
		self.end = end
		self.length = end-start
		self.type=type

def remove_borders(domainsTree,borders):
	print "Removing ",len(borders),"borders"
	for b in borders:
		domainsTree.remove(b)
		
def getDomains(domainsTree,fname,type):
	domains = np.genfromtxt(fname,dtype=np.dtype([('chr','S15'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
	domains = np.sort(domains,order=["chr","start"])
	
	t0 = TAD(domains[0][0],domains[0][1],domains[0][2],type=type)
	for ind,val in enumerate(domains[1:]):
		t = TAD(val[0],val[1],val[2],type=type)
		if t.chr == t0.chr:
			domainsTree[val[1]:val[1]+1] = {"l":t0,"r":t}
		t0 = t

def check_tree(domainsTree):
	for border in domainsTree:
		assert border.data["l"].chr == border.data["r"].chr
		assert border.data["l"].start < border.data["r"].start
		assert border.data["l"].end <= border.data["r"].start
		assert border.data["l"].start < border.data["l"].end
		assert border.data["r"].start < border.data["r"].end

def filter_no_interTADs(domainsTree,maxInterTADlength = 2):
	to_remove = []
	for border in domainsTree:
		if border.data["r"].start - border.data["l"].end > maxInterTADlength:
			to_remove.append(border)
	remove_borders(domainsTree,to_remove)

def filter_same_border(domainsTree,types,accuracy = 2):
	to_remove = []
	for border in domainsTree:
		all_borders = domainsTree[border.begin-accuracy:border.end+accuracy]
		all_types = [b.data["l"].type for b in all_borders if b.data["l"].chr == border.data["l"].chr]
		shared=all([(t in all_types) for t in types])
		if not shared:
			to_remove.append(border)
	
	remove_borders(domainsTree,to_remove)

def filter_only_borders_inside_TAD(domainsTree,ofWhomTADs,distance_from_border = 80000):
	to_leave = []
	print ofWhomTADs
	for border in domainsTree:
		if border.data["l"].type in ofWhomTADs:
			if border.data["l"].length > distance_from_border:
				left = border.data["l"].start
				right = border.data["l"].end
				for b in domainsTree[left:right]:
					if min(b.begin-left,right-b.begin)>distance_from_border and b.data["l"].chr == border.data["l"].chr:
						b.data["l"].within = border.data["l"]
						b.data["r"].within = border.data["l"]
						assert b.data["l"].chr == b.data["l"].within.chr
						assert b.data["r"].chr == b.data["r"].within.chr
						to_leave.append(b)
			if border.data["r"].length > distance_from_border:
				left = border.data["r"].start
				right = border.data["r"].end
				for b in domainsTree[left:right]:
					if min(b.begin-left,right-b.begin)>distance_from_border and b.data["r"].chr == border.data["r"].chr:
						b.data["l"].within = border.data["r"]
						b.data["r"].within = border.data["r"]
						assert b.data["l"].chr == b.data["l"].within.chr
						assert b.data["r"].chr == b.data["r"].within.chr
						to_leave.append(b)

	return IntervalTree(list(set(to_leave)))
	# to_remove = [b for b in domainsTree if not (b in to_leave)]
	# remove_borders(to_remove)

	
def filter_chrm(domainsTree,chromosomes):
	to_remove = []
	for border in domainsTree:
		if not border.data["r"].chr in chromosomes:
			to_remove.append(border)
	remove_borders(domainsTree,to_remove)

def filter_sizes(domainsTree,minsize1=700000,minsize2=1000000,accuracy = 2):
	to_remove = []
	
	if minsize1 > minsize2:
		t=minsize2
		minsize2 = minsize1
		minsize1 = t

	for border in domainsTree:
		all_borders = domainsTree[border.begin-accuracy:border.end+accuracy]
		all_borders = [b for b in all_borders if b.data["l"].chr==border.data["l"].chr]
		pairs = [(min(b.data["r"].length,b.data["l"].length),max(b.data["r"].length,b.data["l"].length)) for b in all_borders]
		if not any([p[0]>minsize1 and p[1]>minsize2 for p in pairs]):
			to_remove.append(border)
			
	remove_borders(domainsTree,to_remove)
	
def filter_gaps(domainsTree):
	to_remove = []
	for border in domainsTree:
		left_contig = converter.chrmTOcontig(border.data["l"].chr,
											border.data["l"].start)
		right_contig = converter.chrmTOcontig(border.data["r"].chr,
											border.data["r"].end)
		if  (left_contig == None) or \
			(right_contig == None) or \
			(left_contig.chrm != right_contig.chrm):
				to_remove.append(border)
			
	remove_borders(domainsTree,to_remove)

def prepE1values(E1_dict):
	return dict([(i,np.genfromtxt(E1_dict[i],dtype=None)) for i in E1_dict])
	
def write_borders(domainsTree,filename,E1):
	domains = list(domainsTree)
	domains = sorted(domains,key=lambda x: (x.data["l"].chr,x.begin))
	with open(filename,"w") as out:
		out.write("\t".join(["Border_chr","Start","End","DomainType",\
							"Within_chr","Start","End","Type",\
							"Left_start","Left_end","Left_E1",\
							"Right_start","Right_end","Right_E1"])+"\n")
		for b in domains:
			E1_key = b.data["l"].type.split("_")[0]
			temp = E1[E1_key]
			E1_l = temp["f2"][(temp["f0"]==b.data["l"].chr) & \
								(temp["f1"] >= b.data["l"].start) & \
								(temp["f1"] <= b.data["l"].end)]
			E1_r = temp["f2"][(temp["f0"]==b.data["l"].chr) & \
							(temp["f1"] >= b.data["r"].start) & \
							(temp["f1"] <= b.data["r"].end)]

			out.write("\t".join(map(str,[\
			b.data["l"].chr,b.begin,b.end,b.data["l"].type,
			b.data["l"].within.chr,b.data["l"].within.start,b.data["l"].within.end,b.data["l"].within.type,
							b.data["l"].start,b.data["l"].end,np.average(E1_l),
							b.data["r"].start,b.data["r"].end,np.average(E1_r)]))+"\n")

base_folder = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/all_domains/"

Domains = {
"ChEF":
	{"Dix":"Dixon-ChEF_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation",
	"Arm":"Arm-ChEF-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation"},#TADtree:""},
"ChME":
	{"Dix":"Dixon-Blood_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation",
	"Arm":"Arm-Blood-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation"},#TADtree""},
"ChIE":
	{"Dix":"Dixon-ChIE_all_HindIII_40k.hm.IC_domains_40KB.jucebox_domains.annotation",
	"Arm":"Arm-ChIE-all-HindIII-40k.hm.gzipped_matrix.jucebox_domains.annotation"},#TADtree""},
}

E1 = {
"ChEF":"/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-100k.hm.eig",
"ChME":"/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-100k.hm.eig",
"ChIE":"/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChIE-all-HindIII-100k.hm.eig"
}
E1=prepE1values(E1)

domainsTree = IntervalTree()
for cell in Domains:
	for alg in Domains[cell]:
		getDomains(domainsTree,base_folder+Domains[cell][alg],cell+"_"+alg)

print len(domainsTree)
#check_tree(domainsTree)
#print len(domainsTree)
print "filter_chrm"
filter_chrm(domainsTree,["chr1","chr2","chr3","chr4","chr11","chr12","chr13","chr14","chr15"])
print len(domainsTree)

print "filter_only_borders_inside_TAD"
within = "ChEF"
distance_from_border = 300000
domainsTree=filter_only_borders_inside_TAD(domainsTree,[within+"_"+i for i in Domains[within].keys()],distance_from_border = distance_from_border)
print len(domainsTree)

print "filter_no_interTADs"
filter_no_interTADs(domainsTree)
print len(domainsTree)

print "filter_gaps"
filter_gaps(domainsTree)
print len(domainsTree)

print "filter_same_border"
same = "ChME"
filter_same_border(domainsTree,[same+"_"+i for i in Domains[same].keys()])
print len(domainsTree)
print "filter_sizes"
filter_sizes(domainsTree)
print len(domainsTree)
write_borders(domainsTree,figure_path+same+"_8eq"+str(distance_from_border/1000)+"kb.txt",E1)