gff_fname = "/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_genomic.gff"

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

		# if geneName == "HBE1" or geneName=="HBE" or geneName=="HBG2"  or geneName=="HBG1":
			# print "\t".join(line)

		
		strand = line[7]
		start = line[3] if contig_type == "chrm" else converter.contig2chrm(line[0],int(line[3]),strand)
		end = line[4] if contig_type == "chrm" else converter.contig2chrm(line[0],int(line[4]),strand)
		if contig_type == "chrm":
			if out_type == "BED":
				out.write(converter.acc2chrm(line[0])+"\t"+start+"\t"+end+"\t"+geneName+"\n")
			elif out_type == "GFF":
				out.write(converter.acc2chrm(line[0])+"\t"+"\t".join(line[1:3])+"\t"+start+"\t"+end+"\t"+"\t".join(line[5:])+"\n")
		elif contig_type=="contig":
			print line
			raise
			# if out_type == "BED":
				# out.write(converter.acc2chrm(line[0])+"\t"+start+"\t"+end+"\t"+geneName+"\n")
			# elif out_type == "GFF":
				# out = open(converter.acc2chrm("\t".join(line[0:3])+start+"\t"+end+"\t"+"\t".join(line[5:]))
			# out.write(start.chrm+"\t"+str(start.nt)+"\t"+str(end.nt)+"\t"+geneName+"\n")
