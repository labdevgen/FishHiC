import sys

def print_usage():
	print "---------------------------"
	print "Usage:"
	print sys.argv[0]," GFF  - for GFF output"
	print sys.argv[0]," BED  - for BED output"
	print sys.argv[0]," FPKM [top|low] path/to/name.fpkm_tracking path/to/out_file  - for FPKM output"
	print sys.argv[0],"for FPKM output format: chrm start end"
	print sys.argv[0],"for FPKM output format: start and end is a region around TSS of top/low25% expressed genes"

if len(sys.argv) == 1:
	print_usage()
	sys.exit()
elif len(sys.argv) == 2:
	if sys.argv[1] == "GFF":
		out_type = "GFF"
	elif sys.argv[1] == "BED":
		out_type = "BED"
	else:
		print_usage()
		sys.exit()
elif len(sys.argv) == 5:
	if sys.argv[1] != "FPKM":
		print_usage()
		sys.exit()
	out_type = "FPKM"
	FPKM_file = sys.argv[3]
	print "reading FPKMs file"
	import numpy as np
	FPKMsArray = np.genfromtxt(FPKM_file,delimiter="\t",dtype=None,skip_header=1)
	gene_names = FPKMsArray["f4"]
	fpkms = FPKMsArray["f9"]
	FPKMs = {}
	for key,val in zip(gene_names,fpkms):
		FPKMs[key] = float(val)
	del FPKMsArray
	
	top = np.percentile(fpkms,95)
	expression_filter_top = [key for key,val in FPKMs.iteritems() if val > top]

	low = np.percentile(fpkms,5)
	if low == 0:
		expression_filter_low = [key for key,val in FPKMs.iteritems() if val <= low]
	else:
		expression_filter_low = [key for key,val in FPKMs.iteritems() if val < low]
	
	if sys.argv[2].lower() == "top":
		expression_filter = expression_filter_top
	elif sys.argv[2].lower() == "low":
		expression_filter = expression_filter_low
		print len(expression_filter)," genes passed FPKM filter"
	else:
		print_usage()
		sys.exit()
else:
	print_usage()
	sys.exit()


from lib_coordinates_converter import *
converter = Ccoordinates_converter(agp_folder = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/AGP",
									chrm2accFile = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
converter.prepare_acc2chrmDict()							
converter.create_agp_dict()							
gff_fname = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff"
if out_type == "BED":
	out = open(gff_fname+".genes.bed","w")
elif out_type == "GFF":
	out = open(gff_fname+".converted.gff","w")
	out.write("#Original GFF file = "+gff_fname+"\n")
	out.write("#Produced by script "+sys.argv[0]+"\n")
elif out_type == "FPKM":
	out = open(sys.argv[4],"w")

not_found_contigs = [] #DEBUG
with open(gff_fname) as gff_file:
	count_noname = 0
	count = 0
	count_notassmebledcontig = 0
	for line in gff_file:
		if line[0] == "#":
			continue
		line = line.split("\t")
		if line[2] != "gene":
			continue
		count += 1
		contig_type = "na"
		if line[0] in converter.acc2chrmDict.keys():
			contig_type = "chrm"
		elif line[0] in converter.contigsInfo.keys():
			contig_type = "contig"
		else:
			not_found_contigs.append(line[0]) #DEBUG
			count_notassmebledcontig += 1
			continue
			
		if len(line) < 9:
			count_noname += 1
			continue
		geneName = line[8].split(";")
		geneName = [i for i in geneName if (i[:5] == "Name=")]
		if len(geneName) != 1:
			count_noname += 1
			print "\t".join(line) #DEBUG
			continue
		geneName = geneName[0][5:]
		
		if out_type=="FPKM" and not (geneName in expression_filter):
			continue
		

		# if geneName == "HBE1" or geneName=="HBE" or geneName=="HBG2"  or geneName=="HBG1":
			# print "\t".join(line)

		
		strand = line[6]
		assert strand in ["+","-"]
		
		
		if out_type=="FPKM":
#			start = line[3] if contig_type == "contig" else converter.chrmTOcontig(line[0],int(line[3]),strand)
#			end = line[4] if contig_type == "contig" else converter.chrmTOcontig(line[0],int(line[4]),strand)
			start = line[3] if contig_type == "chrm" else converter.contig2chrm(line[0],int(line[3]),strand).nt
			start = ChrmPosition(converter.acc2chrm(line[0]),start,strand)
			end = line[4] if contig_type == "chrm" else converter.contig2chrm(line[0],int(line[4]),strand).nt
			end = ChrmPosition(converter.acc2chrm(line[0]),end,strand)
		else:
			start = line[3] if contig_type == "chrm" else converter.contig2chrm(line[0],int(line[3]),strand).nt
			end = line[4] if contig_type == "chrm" else converter.contig2chrm(line[0],int(line[4]),strand).nt
		
		if start.nt > end.nt:
			if converter.contigsInfo[start.chrm][3] != "-":
				print start
				print end
				print "\t".join(line)
				print converter.contigsInfo[start.chrm]
				raise
			start0 = start
			start = end
			end = start0
		
		if strand == "+":
			end.nt = start.nt+5000
			start.nt = max(0,start.nt-5000)
		else:
			end.nt = end.nt+5000
			start.nt = max(0,end.nt-5000)
		
		if contig_type == "chrm":
			if out_type == "BED":
				out.write(converter.acc2chrm(line[0])+"\t"+start+"\t"+end+"\t"+geneName+"\n")
			elif out_type == "GFF":
				out.write(converter.acc2chrm(line[0])+"\t"+"\t".join(line[1:3])+"\t"+start+"\t"+end+"\t"+"\t".join(line[5:])+"\n")
			elif out_type == "FPKM":
				if start.nt == 5156723:
					print "\t".join(line)
				out.write(start.chrm+"\t"+str(start.nt)+"\t"+str(end.nt)+"\n")
		elif contig_type=="contig":
			print line
			raise
			# if out_type == "BED":
				# out.write(converter.acc2chrm(line[0])+"\t"+start+"\t"+end+"\t"+geneName+"\n")
			# elif out_type == "GFF":
				# out = open(converter.acc2chrm("\t".join(line[0:3])+start+"\t"+end+"\t"+"\t".join(line[5:]))
			# out.write(start.chrm+"\t"+str(start.nt)+"\t"+str(end.nt)+"\t"+geneName+"\n")

out.close()
print count_noname," genes scipped (no gene name)"
print float(count_notassmebledcontig*100)/count,"% of genes scipped (located in not assmebled contigs)"