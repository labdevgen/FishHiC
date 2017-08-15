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
##############################

chrms_axes=False # In this var stored information whether chromosomes are aready drawn

bincounts=[5]+[20]+[40]+[60]+[80]+[100]
color_count=0
for bincount in bincounts:
	print "----bincount=",bincount,"--------"
#w=raw_input("Enter file names, n=next ")   #Uncoment this line and
#	w="n" 										#remove this line to allow user to enter filenames interactively

	filenames={}
	dis={}

#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_Fib_full.hdf5.matrix.di']='FibDI'
#filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_Sp_full.hdf5.matrix.di']='SpDI'
#filenames['/mnt/storage/home/vsfishman/HiC/data/raw/heatmap-res-100KB_Fib_full.hdf5.matrix.di']='Fib-rawDI'
#filenames['/mnt/storage/home/vsfishman/HiC/data/raw/heatmap-res-100KB_Sp_full.hdf5.matrix.di']='Sp-rawDI'

#	filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_Fib_full.hdf5']='FibMatr'
#	filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_Sp_full.hdf5']='SpMatr'
	filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_ESC_full.hdf5']='EscFullMatr'
	filenames['/mnt/storage/home/vsfishman/HiC/data/heatmap-res-100KB_SRR443884_mESC_2.hdf5']='Esc2Matr'

#filenames['/mnt/storage/home/vsfishman/HiC/data/raw/heatmap-res-100KB_Fib_full.hdf5.matrix']='Fib-rawMatr'
#filenames['/mnt/storage/home/vsfishman/HiC/data/raw/heatmap-res-100KB_Sp_full.hdf5.matrix']='Sp-rawMatr'



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
					genome_db = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])
					zcount=0 # count of zerro lines in dataset
				BD = binnedData.binnedData(resolution, genome_db)
				BD.simpleLoad(i, 'HindIII_GM_1')
				for j in range(len(BD.dataDict['HindIII_GM_1'])):
					cur_chr_indx=genome_db.chrmIdxBinCont[j]
		#			print "for bin ",j," appending string ",genome_db.chrmStartsBinCont[cur_chr_indx]," to ",genome_db.chrmEndsBinCont[cur_chr_indx]
					start0=genome_db.chrmStartsBinCont[cur_chr_indx]
					end0=genome_db.chrmEndsBinCont[cur_chr_indx]
		#			start=start0+(end0-start0)/4
		#			end=end0-(end0-start0)/4
					#start=0
					#end=len(BD.dataDict['HindIII_GM_1'][j])-1
			
					start=max(start0,j-bincount)
					end=min(end0,j+bincount)
					dis[i].append(BD.dataDict['HindIII_GM_1'][j,start:end].tolist())

				if len(chrms_in_bins)==0:
					for j in range(len(BD.dataDict['HindIII_GM_1'])):
						chrms_in_bins.append(genome_db.chrmIdxBinCont[j])
				del BD
	

	w=50
	list_analize_type="Sp"
	
	filenames_keylist=sorted(filenames.keys())
	
	while (w != "e") and (w != "E"):
		compare_list={}
		
		inner_window_size=None

		for i in range(0,len(filenames_keylist)):
				for j in range (i+1,len(filenames_keylist)):
					print filenames[filenames_keylist[i]]+"_vs_"+filenames[filenames_keylist[j]]+"="+str(i*j+j)
					compare_list[str(i*j+j)]=[]
					compare_list[str(i*j+j)].append(filenames_keylist[i])
					compare_list[str(i*j+j)].append(filenames_keylist[j])
			

		q="1"
		filename_addition = ""
		for i in q.strip().split():
			write_cor_to_file(dis[compare_list[i][0]],dis[compare_list[i][1]],
									o_filename,chrms_in_bins,chrms_axes,w,False,
									filenames[compare_list[i][0]]+" vs "+filenames[compare_list[i][1]]+"_"+str(bincount),
									color_count,list_analize_type)
			filename_addition += filenames[compare_list[i][0]]+"_vs_"+filenames[compare_list[i][1]]
			if (color_count >= len(colors)-1):
				color_count=0
			else:
				color_count+=1
			chrms_axes = True

		w="e"
plt.subplots_adjust(bottom=0.15)
o_filename_complete="Overviw-ESCvsESC.png"
f=open(o_filename_complete,"wb")
plt.savefig(o_filename_complete,format='png',dpi=600,bbox_inches='tight')	
f.close()
plt.clf()
