from intervaltree import Interval, IntervalTree
import numpy as np

from mirnylib.systemutils import setExceptionHook
setExceptionHook()


sp1 = "hs"
sp1domainsFile = "hs.total.combined.domain"
sp2 = "mm"
sp2domainsFile = "mm.total.HindIII.combined.domain"

orthologsFile = "HumanGRCh38.p7-mouseGRCm38.p4-orthologs-mart_export.txt"
syntFile = "mouse_human.synteny.txt"
#syntFile = "test_mm_hs.txt"

def get_synteny(fname):
	synteny = np.genfromtxt(fname,dtype=np.dtype(
									[('sp2name','S50'),
									('sp2region','S100'),
									('sp1name','S50'),
									('sp1region','S100')
									]))
	byChr_Synteny = {}
	byChr_Synteny[sp1] = {}
	for region in synteny:
		sp1chrm,sp1start,sp1end = tuple(region["sp1region"].split(":")[2:5])
		sp1chrm = "chr"+sp1chrm
		sp1start,sp1end = int(sp1start),int(sp1end)
		
		sp2chrm,sp2start,sp2end =  tuple(region["sp2region"].split(":")[2:5])
		sp2chrm = "chr"+sp1chrm
		sp2tart,sp1end = int(sp1start),int(sp1end)

		if not sp1chrm in byChr_Synteny[sp1].keys():
			byChr_Synteny[sp1][sp1chrm] = IntervalTree()
		sp1Interval = Interval(sp1start,sp1end,
								{"chr":sp2chrm,"start":sp2start,"end":sp2end})
		byChr_Synteny[sp1][sp1chrm].add(sp1Interval)
	return byChr_Synteny

def get_domains(fname):
	domains = np.genfromtxt(fname,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
	domains_dict = {}
	print "Generating intervals object"
	count = 0
	for domain in domains:
		count += 1 
		if not domain["chrm"] in domains_dict.keys():
			domains_dict[domain["chrm"]]=IntervalTree()
			prev_domain = None
		if prev_domain==None:
			domains_dict[domain["chrm"]].add(Interval(domain["start"],domain["end"],count))
		else:
			domains_dict[domain["chrm"]].add(Interval(prev_domain["start"],domain["end"],count))
	
		prev_domain = domain
	
	return domains_dict

def get_domains2(fname):
	domains = np.genfromtxt(fname,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
	domains = np.sort(domains,order=["chrm","start"])
	
	domains_dict = {}
	print "Generating intervals object"
	count = 0
	for domain in domains:
		count += 1 
		if not domain["chrm"] in domains_dict.keys():
			domains_dict[domain["chrm"]]=IntervalTree()
			prev_domain = None
		if prev_domain!=None:
			count += 1
			domains_dict[domain["chrm"]].add(Interval(prev_domain["start"]+(prev_domain["end"]-prev_domain["start"])/2,
											domain["start"]+(domain["end"]-domain["start"])/2,count))
		prev_domain = domain
	
	return domains_dict

	
def get_boundaries(fname):
	domains = np.genfromtxt(fname,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
	domains_dict = {}
	print "Generating intervals object"
	count = 0
	for domain in domains:
		count += 1 
		if not domain["chrm"] in domains_dict.keys():
			domains_dict[domain["chrm"]]=IntervalTree()
			prev_domain = None
		domains_dict[domain["chrm"]].add(Interval(domain["start"],domain["end"],count))
		if prev_domain != None:
			domains_dict[domain["chrm"]].add(Interval(prev_domain["end"],domain["start"],{"ID":count,"Type":"border"}))
		else:
			domains_dict[domain["chrm"]].add(Interval(prev_domain["start"],domain["end"],{"ID":count,"Type":"domain"}))
		prev_domain = domain
	
	return domains_dict

	
def generate_orthologs(fname):
	# Ensembl Gene ID	++++++++++++0
	# Ensembl Transcript ID	1
	# Mouse Ensembl Gene ID	2
	# Mouse associated gene name	3
	# Mouse Chromosome Name	4
	# Mouse Chromosome start (bp)	5
	# Mouse Chromosome end (bp)	6
	# Mouse homology type	7
	# Chromosome Name	8
	# Gene Start (bp)	9
	# Gene End (bp)	10

	orthologs = np.genfromtxt(fname,dtype=np.dtype(
									[('sp1geneID','S50'),
									('sp2geneName','S50'),
									('sp2chrm','S50'),
									('sp2start',np.uint32),
									('sp2end',np.uint32),
									('homolType','S50'),
									('sp1chrm','S50'),
									('sp1start',np.uint32),
									('sp1end',np.uint32)]),usecols = (0,3,4,5,6,7,8,9,10)) 
									
	orthologs = orthologs[orthologs["homolType"]=="ortholog_one2one"]
	
	byChr_orthologs={}
	byChr_orthologs[sp1]={}
	byChr_orthologs[sp2]={}
	gene_IDs = []
	for gene in orthologs:
		if gene["sp1geneID"] in gene_IDs:
			continue
		if min(gene["sp1end"]-gene["sp1start"],gene["sp2end"]-gene["sp2start"])<5000:
			continue
		sp1chrm = "chr"+gene['sp1chrm']
		sp2chrm = "chr"+gene['sp2chrm']
		if not sp1chrm in byChr_orthologs[sp1].keys():
			byChr_orthologs[sp1][sp1chrm] = IntervalTree()
		sp1Interval = Interval(gene["sp1start"],gene["sp1end"],
								{"name":gene["sp1geneID"],"or_chr":sp2chrm,"or_start":gene["sp2start"],"or_end":gene["sp2end"],"or_name":gene["sp2geneName"]})
		byChr_orthologs[sp1][sp1chrm].add(sp1Interval)
	return byChr_orthologs

	
print "Reading domains1"
sp1Domains = get_domains2(sp1domainsFile)
print "Reading domains2"
sp2Domains = get_domains2(sp2domainsFile)
print "Readsing Syntenic regions"
synt_regions = get_synteny(syntFile)
print "Reading orthologs"
orthologs = generate_orthologs(orthologsFile)


plus_plus = 0
plus_minus = 0
total_genes = 0
for sp1chr in synt_regions[sp1]:
	for sp1SyntRegion in synt_regions[sp1][sp1chr]:
		print "Region: ",sp1chr,":",sp1SyntRegion.begin,"-",sp1SyntRegion.end
		if not sp1chr in orthologs[sp1].keys(): #no genes on this chrm
			print sp1chr,orthologs[sp1].keys()
			continue
		sp1genes = list(orthologs[sp1][sp1chr][sp1SyntRegion.begin:sp1SyntRegion.end])
		print "Region contains ",len(sp1genes)," genes"
#		matrix = np.zeros(shape=(len(sp1genes),len(sp1genes)))
		total_genes += len(sp1genes)
		for i in xrange(len(sp1genes)):
			for j in xrange(i+1,len(sp1genes)):
				gene1 = sp1genes[i]
				gene2 = sp1genes[j]
				domains1 = [q.data for q in sp1Domains[sp1chr][gene1.begin:gene1.end]]
				domains2 = [q.data for q in sp1Domains[sp1chr][gene2.begin:gene2.end]]
				sp1overlap = len(set(domains1+domains2))<(len(domains1)+len(domains2))
				if sp1overlap:
					# print "Overlap in human:",domains1,domains2
					# for d in sp1Domains[sp1chr]:
						# if d.data in domains1+domains2:
							# print "--",d 
						
					sp2chr = gene1.data["or_chr"]
					if sp2chr!=gene2.data["or_chr"]:
						continue
					domains1 = [q.data for q in sp2Domains[sp2chr][gene1.data["or_start"]:gene1.data["or_end"]]]
					domains2 = [q.data for q in sp2Domains[sp2chr][gene2.data["or_start"]:gene2.data["or_end"]]]
					sp2overlap = len(set(domains1+domains2))<(len(domains1)+len(domains2))
					if sp1overlap and sp2overlap:
						plus_plus += 1
					else:
						plus_minus += 1
						# print "-------------"
						# print gene1,gene2
						# print "-------------"
	print plus_plus,plus_minus
print total_genes,plus_plus,plus_minus,plus_plus/float(plus_plus+plus_minus)
