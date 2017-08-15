from lib_06_DI_Sperman import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("o_file_name")
args = parser.parse_args()

o_filename=args.o_file_name

chrms_in_bins=[] #Item number i of this array is chromosome number for bin number i
genome_name='mm9'
genome_db=None
resolution=0

##############################
list_analize_type="proc"
#possible values: "Sp" - for Spearman cr
#					"Evkl" - for Evklid metrik
#					"proc" - percentage of bins, with wich they have same interactions
#					"pears" - Pearson correlation				
#		Add _w to analizis_type to preform it with inner window, e.g. Sp_w is Spearman cr with inner window
##############################

#w=raw_input("Enter file names, n=next ")   #Uncoment this line and
w="n" 										#remove this line to allow user to enter filenames

filenames={}
dis={}

#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_Fib_full.hdf5.matrix.di']='FibDI'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_Sp_full.hdf5.matrix.di']='SpDI'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_ESC_full.hdf5.matrix.di']='ESC_full'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_SRR443884_mESC_2.hdf5.matrix.di']='ESC'

filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_Fib_full2.hdf5']='Fibroblasts'
filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_Sp_full2.hdf5']='Sperm cells'
filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_ESC_full.hdf5']='ESC'
filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-1000KB_SRR443884_mESC_2.hdf5']='ESC_2'

#filenames['/mnt/storage/home/vsfishman/HiC/data/raw/heatmap-res-100KB_Fib_full.hdf5.matrix']='Fib-rawMatr'
#filenames['/mnt/storage/home/vsfishman/HiC/data/raw/heatmap-res-100KB_Sp_full.hdf5.matrix']='Sp-rawMatr'


while (w != "n") and (w != "N"):
  filenames[w]=''
  q=raw_input("Enter name for file "+w+' ')
  filenames[w]=str(q).strip()
  w=raw_input("Enter file names, n=next ")   

if len(filenames.keys())<2:
   print "Please provide at least 2 files"
   exit()
   
for i in filenames.keys():  
	dis[i]=[]
	print "Reading file "+i+"\n"

	if (i.split('.')[-1]=='hdf5'):
		if (resolution==0):
			raw_heatmap = h5dict.h5dict(i, mode='r')
			resolution = int(raw_heatmap['resolution'])
			del raw_heatmap
		if (genome_db==None):
			genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
		zcount=0 # count of zerro lines in dataset
		BD = binnedData.binnedData(resolution, genome_db)
		BD.simpleLoad(i, 'HindIII_GM_1')
		scan_genome_from_bin=0
#		scan_genome_from_bin=genome_db.chrmStartsBinCont[4]
		for j in range(scan_genome_from_bin,len(BD.dataDict['HindIII_GM_1'])):
#		for j in range(scan_genome_from_bin,genome_db.chrmEndsBinCont[0]):
			cur_chr_indx=genome_db.chrmIdxBinCont[j]
#			if cur_chr_indx > 0: continue
#			print "for bin ",j," appending string ",genome_db.chrmStartsBinCont[cur_chr_indx]," to ",genome_db.chrmEndsBinCont[cur_chr_indx]
			start0=genome_db.chrmStartsBinCont[cur_chr_indx] #start is the begining of current chromosome
			end0=genome_db.chrmEndsBinCont[cur_chr_indx] #start is the end of current chromosome
#			start=start0+(end0-start0)/4
#			end=end0-(end0-start0)/4
			start=max(j-90,start0) #before this constatn was 50
			end=min(j+90,end0)
#			end=len(BD.dataDict['HindIII_GM_1'][j])-1
#			bincount=10
#			start=max(start0,j-bincount)
#			end=min(end0,j+bincount)
			dis[i].append(BD.dataDict['HindIII_GM_1'][j,start0:end0].tolist())
			s=sum(dis[i][j-scan_genome_from_bin])
			if (s==0):
			  zcount += 1
		print "Number of zerro lines=",zcount
		if len(chrms_in_bins)==0:
			for j in range(len(BD.dataDict['HindIII_GM_1'])):
				chrms_in_bins.append(genome_db.chrmIdxBinCont[j])
		print len(chrms_in_bins)
		del BD
	elif (i.split('.')[-1]=='matrix'):
			f=open(i,'r')
			for line in f:
				s=[]
				splited_line=line.split("\t")
				for t in splited_line[3:]:
					s.append(float(t))
				dis[i].append(s)
				chrms_in_bins.append(splited_line[0])
			f.close()
		# Variant of analizis to correlate total sums
#		import os
#		if os.path.isfile(i+'-sum'):
#			print "ustig file "+i+'-sum'
#			for line in f:
#				splited_line=line.split("\t")
#				dis[i].append(splited_line[3])
#				chrms_in_bins.append(splited_line[0])
#		else:
#			q=raw_input("Write output data to file? (y/n) ")
#			WriteOut=False
#			if ((q=='y') or (q=='Y')):
#				f_out=open(i+'-sum','w')
#				WriteOut=True
#			for line in f:
#				s=0
#				splited_line=line.split("\t")
#				for t in splited_line[3:]:
#					s+=float(t)
#				dis[i].append(s)
#				print s 
#				chrms_in_bins.append(splited_line[0])
#				if WriteOut:
#					f_out.write(splited_line[0]+"\t"+splited_line[1]+"\t"+splited_line[2]+"\t"+str(s)
#			if WriteOut:
#				f_out.close()
	elif (i.split('.')[-1]=='di'):			
		f=open(i,'r')
		for line in f:
			dis[i].append(float(line.split()[3]))
			chrms_in_bins.append(line.split()[0])
		f.close()
	print "Length of data in file is "+str(len(dis[i]))
	


chrms_axes=False # In this var stored information whether chromosomes are aready drawn

filenames_keylist=sorted(filenames.keys())

analysis_types = ["Sp","Evkl","pears"]
window_sizes = [1,1,1]

#w=raw_input("Enter window size, e=exit ")
#list_analize_type=raw_input("Enter list_analize_type ").strip()
#while (w != "e") and (w != "E"):

for interation_index in [0,1,2]:
  w = window_sizes[interation_index]
  list_analize_type = analysis_types[interation_index]
  
  compare_list={}
  color_count=0
  inner_window_size=None

  for i in range(0,len(filenames_keylist)):
		for j in range (i+1,len(filenames_keylist)):
			print filenames[filenames_keylist[i]]+"_vs_"+filenames[filenames_keylist[j]]+"="+str(i*j+j)
			compare_list[str(i*j+j)]=[]
			compare_list[str(i*j+j)].append(filenames_keylist[i])
			compare_list[str(i*j+j)].append(filenames_keylist[j])
			

  q='2 1 3 6'
#  q='6'
#  q=raw_input("Which graphs would you like to show, a=all ")
  filename_addition = ""
  if (q.strip()=='a') or (q.strip()=='A'):
	for i in range(0,len(filenames_keylist)):
		print "Plotting picture ", i+1, "( of ",len(filenames_keylist),")"
		for j in range (i+1,len(filenames_keylist)):
			write_cor_to_file(dis[filenames_keylist[i]],dis[filenames_keylist[j]],
									o_filename,chrms_in_bins,chrms_axes,genome_db,w,False,
									filenames[filenames_keylist[i]]+" vs "+filenames[filenames_keylist[j]],
									color_count,list_analize_type)
			chrms_axes = True
			filename_addition += filenames[filenames_keylist[i]]+"_vs_"+filenames[filenames_keylist[j]]
  else:
	do_correction='y'
#	do_correction=raw_input("Do correction(y/n)? ")
	correction_list=[]
	if do_correction == "y":
		do_correction = True
		cor_ind='2'
#		cor_ind=raw_input("Enter correction set number ")
		q1=""
		for a in q.strip().split():
			if a != cor_ind.strip().split()[0]: #remove correction set from list
				q1=q1+a+" " 
		q=q1
	else:
		do_correction  = False
		cor_ind=''
		
	for ind,i in enumerate(cor_ind.strip().split()+q.strip().split()):
			print "Plotting picture ", i
			if (do_correction and ind==0):
				correction_list=write_cor_to_file(dis[compare_list[i][0]],dis[compare_list[i][1]],
									o_filename,chrms_in_bins,chrms_axes,genome_db,w,False,
									filenames[compare_list[i][0]]+" vs "+filenames[compare_list[i][1]],
									color_count,list_analize_type,False,correction_list)
			else:
				write_cor_to_file(dis[compare_list[i][0]],dis[compare_list[i][1]],
									o_filename,chrms_in_bins,chrms_axes,genome_db,w,False,
									filenames[compare_list[i][0]]+" vs "+filenames[compare_list[i][1]],
									color_count,list_analize_type,do_correction,correction_list)
									
			filename_addition += filenames[compare_list[i][0]]+"_vs_"+filenames[compare_list[i][1]]
			if (color_count >= len(colors)-1):
				color_count=0
			else:
				color_count+=1
			chrms_axes = True
  o_filename_complete=o_filename+str(w)+filename_addition+"_1MB"+"."+list_analize_type+'.png'
  if do_correction:
	  o_filename_complete += '.corrected'
  o_filename_complete += '.png'
  
  plt.subplots_adjust(bottom=0.15)
  f=open(o_filename_complete,"wb")
  plt.savefig(o_filename_complete,format='png',dpi=800,bbox_inches='tight')	
  f.close()
  plt.clf()
  chrms_axes = False
#  w=raw_input("Enter window size, e=exit ")
#  list_analize_type=raw_input("Enter list_analize_type ").strip()