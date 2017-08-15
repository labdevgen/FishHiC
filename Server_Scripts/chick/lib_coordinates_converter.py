import os

#This module was designed to convert fasta names between different formats
#formats are:
#chrm examples "chr1","chrZ"
#acc example "NC_006088.4"
#contig example "NT_455768.1"
#To convert different formats we need to types of information
#first is AGP files describing relationships between assembled chromosomes and contigsInfo
#AGP files suppoused to be in AGP directory and have extansion .agp
#AGP file example:
	##comment lines
	#gap_type      linkage        evidence
	#NC_006097.4     1       500000  1       N       500000  centromere      no     na
	#NC_006097.4     500001  632397  2       O       NT_455918.1     1       132397 +
	#NC_006097.4     632398  632497  3       U       100     contig  no      na
	#NC_006097.4     632498  678831  4       O       NT_455919.1     1       46334  +
#second information is chrm2accFile
#example:
##Chromosome     Accession.version
#1       NC_006088.4
#2       NC_006089.4
#3       NC_006090.4
#It is recommended to create new instance of class Ccoordinates_converter and imideatly supply both AGP and chrm info
#		converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
#									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
#		converter.prepare_acc2chrmDict()
#		converter.create_agp_dict()



class ChrmPosition:
	def __init__(self,chrm,nt,strand):
		self.chrm = chrm
		self.nt = int(nt)
		self.strand = strand
	def __str__(self):
		return "\t".join(map(str,[self.chrm,self.nt,self.strand]))
	
class Ccoordinates_converter:
	def getChrmsByType(self,format="chrm"):
	#input: 
	#format - string, one of "chrm","acc","contig"
	#chrm - return chrms like "chr1","chrZ"
	#acc -  return chrms like "NC_006088"
	#contig - return contigs like "NT_455768"
	#output:
	#dictionary with keys "micro","macro","sex"
	#under each key list of names (according to format) with corresponding chromosomes
		assert format in ["chrm","acc","contig"]
		if format == "chrm":
			return {"micro":self._MicroChrms,"macro":self._MacroChrms,"sex":self._SexChrms}
		if format == "acc":
			if self.acc2chrmDict==None:
				self.prepare_acc2chrmDict()
			return {"micro":[self.chrm2acc(i) for i in self._MicroChrms],
					"macro":[self.chrm2acc(i) for i in self._MacroChrms],
					"sex":[self.chrm2acc(i) for i in self._SexChrms]}
		if format == "contig":
			if self.chrm2contig==None:
				self.create_agp_dict()
			return {"micro":[i for chr in self._MicroChrms for i in self.chrm2contig[chr]],
					"macro":[i for chr in self._MacroChrms for i in self.chrm2contig[chr]],
					"sex":[i for chr in self._SexChrms for i in self.chrm2contig[chr]]}

		
	def __init__(self,agp_folder = None, chrm2accFile = None):
		self.agp_folder = agp_folder
		self.chrm2accFile = chrm2accFile
		self._MacroChrms = ["chr"+str(i) for i in range(1,11)]
		self._MicroChrms = ["chr"+str(i) for i in range(11,34) if not i in [29]]+["chrLGE64"]
		self._SexChrms = ["chrZ","chrW"]

		
	def setAGPfolder(self,agp_folder):
		self.agp_folder = agp_folder
		
	def setChrm2accFile(self,chrm2accFile):
		self.chrm2accFile = chrm2accFile

	def prepare_acc2chrmDict(self,chrm2accFile=None):
		if chrm2accFile==None:
			if self.chrm2accFile == None:
				raise Exception("Please provide folder with chromosome->acc info when running this function first time")
		else:
			self.chrm2accFile = chrm2accFile
				
		with open(self.chrm2accFile) as f:
			self.chrm2accDict = dict(("chr"+line.strip().split()[0],line.strip().split()[1]) for line in f if line[0] != "#")
			f.seek(0)
			self.acc2chrmDict = dict((line.strip().split()[1],"chr"+line.strip().split()[0]) for line in f if line[0] != "#")
		
	def acc2chrm(self,acc,chrm2accFile = None):
	#input - acc name of chromosome
	#output - chr name of chromosme
		if self.acc2chrmDict == None:
			self.prepare_acc2chrmDict(chrm2accFile)
		return self.acc2chrmDict[acc]

	def chrm2acc(self,chrm,chrm2accFile = None):
		if self.chrm2accDict == None:
			self.prepare_acc2chrmDict(chrm2accFile)
		return self.chrm2accDict[chrm]
		
		
	def create_agp_dict(self,agp_folder = None):
		if agp_folder == None:
			if self.agp_folder == None:
				raise Exception("Please provide folder with agp files")
			else:
				agp_folder = self.agp_folder
		else:
			self.agp_folder = agp_folder
			
		print "Scanning AGP folder " + agp_folder
		chrm2contig = {}
		contigsInfo = {}
		chrmsLens = {}
		
		agp_files = [f for f in os.listdir(agp_folder) if f.endswith(".agp")]
		for f in sorted(agp_files,key=len):
			with open(agp_folder+"/"+f) as agp_file:
				chrName = f.split(".")[0]
				chrm2contig[chrName] = []
				chr_len = 0
				for line in agp_file:
					if line[0]=="#": continue
					line = line.split()
					if line[4] == "O" or line[4] == "W": #O is a contig, U is a gap, N is a centromere, W - contig on chrm W
							contigName = line[5]
							assert int(line[1])<int(line[2])
							contigsInfo[contigName] = [chrName,int(line[1]),int(line[2]),line[-1]] #contigName --> [chr,contigStartBp,contigEndBp,orientataion]
							chrm2contig[chrName].append(contigName)
					assert chr_len < int(line[2])
					chr_len = int(line[2])
				chrmsLens[chrName] = int(chr_len)
			
		self.chrm2contig = chrm2contig
		self.contigsInfo = contigsInfo
		self.chrmsLens = chrmsLens
		print "Done"
	
	def contig2chrm(self,contig,nt,strand="na"):
	#input - position on chromosome, chrm name in contig format, nt according to position on contig
	#output - position on chromosome, chrm name in chrm format,  nt according to position on whole chromosome
		assert strand in ["na","+","-"]
		if self.contigsInfo == None:
			self.create_agp_dict()
		
		if not (contig in self.contigsInfo.keys()):
			raise Exception("Contig not found in agpInfo")
		chrm = self.contigsInfo[contig][0]
		if self.contigsInfo[contig][3] == "+":
			nt = self.contigsInfo[contig][1] + nt - 1
		else:
			nt = self.contigsInfo[contig][2] - nt + 1
			if strand != "na":
				strand = "-" if strand=="+" else "-"
			
		assert nt > 0
		assert nt <= self.chrmsLens[chrm]
		
		return ChrmPosition(chrm,nt,strand)
	
	def chrmTOcontig(self,chrm,nt,strand="na"):
	#input - position on chromosome, chrm name in chrm or acc format, nt according to position on whole chromosome
	#output - position on chromosome, chrm name in contig format, nt according to position on defined contig
		assert strand in ["+","-","na"]
		if self.contigsInfo == None:
			self.create_agp_dict()

		if not chrm in self.chrm2contig.keys():
			if self.acc2chrmDict == None:
				self.prepare_acc2chrmDict(chrm2accFile)
			if not chrm in self.acc2chrmDict.keys():
				raise Exception("Chromosome "+str(chrm)+" not found")
			else:
				chrm = self.acc2chrmDict[chrm]
		
		found = False
		for contig in self.chrm2contig[chrm]:
			contigName = contig
			contig = self.contigsInfo[contig]
			if contig[1]<=nt<=contig[2]:
				found = True
				break
		
		if not found:
			print "WARNING: contig not found for coordinate ",chrm,"\t",nt
			return None
			
		assert contig[1]<=nt<=contig[2]
		assert contig[0] == chrm
		
		if contig[3] == "+":
			nt = nt - contig[1] + 1
		elif contig[3] == "-":
			nt = contig[2] - nt + 1
			if strand != "na":
				strand = "-" if strand=="+" else "-"
		assert nt > 0
		assert nt <= self.chrmsLens[chrm]
		return ChrmPosition(contigName,nt,strand)