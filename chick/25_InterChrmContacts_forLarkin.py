import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import sys
import numpy as np

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData

from itertools import cycle

from mirnylib.systemutils import setExceptionHook
setExceptionHook()

genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

heatmap_filepath=sys.argv[1]
resolution = int(heatmap_filepath.split("-")[-1].split("k")[0])*1000
print "Resolution determined: ",resolution

specialTransFile=sys.argv[2]

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
import figPath
import ntpath
figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"_"+'kb_Trans-EBR_regions.png'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad(heatmap_filepath, 'heatmap')
q=BD.dataDict['heatmap']

all_trans = []
for i in range(genome_db.chrmCount-1): # minus 1 is because all contacts of last chr will be counted at the end of for
	all_trans+=list(
					q[genome_db.chrmStartsBinCont[i]:genome_db.chrmEndsBinCont[i],
					  genome_db.chrmStartsBinCont[i+1]:].flatten()
					)


fig1=plt.figure(1)
plt.boxplot(all_trans,whis=[1,99],showfliers=False)
x_labels = ["","All trans (genome)",""]

plt.title(heatmap_filepath+"\n"+specialTransFile)

specialTrans = dict([i.strip().split(":") for i in open(specialTransFile).readlines()])

colors = cycle(["k","b","g","r"])
markers = cycle(["o","v","x","*","8"])
legends=""
to_plot={}
to_plot2={}

for x_ind,t in enumerate(sorted(specialTrans.keys())):
	it = iter(t.split())
	set1 = zip(it,it,it) # split list to values chr-start-end
	it = iter(specialTrans[t].split())
	set2 = zip(it,it,it)

	for reg1 in set1:
		color = colors.next()
		for reg2 in set2:
			chr1,start1,end1 = reg1[0],reg1[1],reg1[2]
			chr2,start2,end2 = reg2[0],reg2[1],reg2[2]
			assert chr1 != chr2
			chr1,chr2 = genome_db.label2idx[chr1],\
						genome_db.label2idx[chr2]
			start1,start2,end1,end2 = max(0,int(start1)/resolution),\
									max(0,int(start2)/resolution),\
									min(genome_db.chrmLensBin[chr1],int(end1)/resolution),\
									min(genome_db.chrmLensBin[chr2],int(end2)/resolution)
			current_trans = q[genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1,
							genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1].flatten()
			to_plot["_".join([genome_db.idx2label[chr1],
								genome_db.idx2label[chr2]])
					] = q[genome_db.chrmStartsBinCont[chr1]:genome_db.chrmStartsBinCont[chr1+1],
												genome_db.chrmStartsBinCont[chr2]:genome_db.chrmStartsBinCont[chr2+1]]
			marker = markers.next()
			
			#Trans for region 1/2
			plt.plot([len(x_labels)]*len(current_trans),current_trans,marker=marker,color=color,ls="None",
				label=" ".join([genome_db.idx2label[chr1],str(start1*resolution),genome_db.idx2label[chr2],str(start2*resolution)])
				)
			x_labels.append(genome_db.idx2label[chr1]+" "+str(start1*resolution/1000)+" vs "+genome_db.idx2label[chr2]+" "+str(start2*resolution/1000))
			
			#All trans for region1
			current_trans_left = q[genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1,
							0:genome_db.chrmStartsBinCont[chr1]-1].flatten()
			current_trans_right = q[genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1,
							genome_db.chrmEndsBinCont[chr1]+1:genome_db.chrmEndsBinCont[-1]].flatten()
			current_trans = np.concatenate((current_trans_left,current_trans_right))
			plt.boxplot(current_trans,positions=[len(x_labels)],whis=[1,99],showfliers=False)
			x_labels.append(genome_db.idx2label[chr1]+"_"+str(start1*resolution/1000)+"k")
			
			#All trans for region2
			current_trans_left = q[genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1,
							0:genome_db.chrmStartsBinCont[chr2]-1].flatten()
			current_trans_right = q[genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1,
							genome_db.chrmEndsBinCont[chr2]+1:genome_db.chrmEndsBinCont[-1]].flatten()
			current_trans = np.concatenate((current_trans_left,current_trans_right))
			plt.boxplot(current_trans,positions=[len(x_labels)],whis=[1,99],showfliers=False)
			x_labels.append(genome_db.idx2label[chr2]+"_"+str(start2*resolution/1000)+"k")
			
			#Plot of all trans of region1
			current_trans = q[genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1,
							:]
			current_trans[:,genome_db.chrmStartsBinCont[chr1]:genome_db.chrmEndsBinCont[chr1]] = 0
			current_trans2 = q[genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1,
							genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1].flatten()
			assert len(current_trans2.flatten())==1
			current_trans2=current_trans2[0]
			fig2=plt.figure(2)
			plt.subplot(211)
			plt.title(genome_db.idx2label[chr1]+"_"+str(start1*resolution/1000),fontsize="xx-small")
			plt.plot(range(len(q)),current_trans[0],label=i,marker=".",color="k",ls="None",markersize=0.4)
			plt.plot(genome_db.chrmStartsBinCont[chr2]+start2,current_trans2,marker="o",color="b",markersize=2)
			plt.gca().set_xticks(genome_db.chrmStartsBinCont)
			plt.gca().set_xticklabels([genome_db.idx2label[i] for i in range(genome_db.chrmCount)],rotation=90,fontsize="xx-small")
#			fig2.savefig(figure_path+"."+genome_db.idx2label[chr1]+"_"+str(start1*resolution/1000)+".png",dpi=700)
			plt.figure(1)

			#Plot of all trans of region2
			current_trans = q[genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1,
							:]
			current_trans[:,genome_db.chrmStartsBinCont[chr2]:genome_db.chrmEndsBinCont[chr2]] = 0

			current_trans2 = q[genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1,
							genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1].flatten()
			assert len(current_trans2.flatten())==1
			current_trans2=current_trans2[0]
			fig2=plt.figure(2)
			plt.subplot(212)
			plt.title(genome_db.idx2label[chr2]+"_"+str(start2*resolution/1000),fontsize="xx-small")
			plt.plot(range(len(q)),current_trans[0],label=i,marker=".",color="k",ls="None",markersize=0.4)
			plt.plot(genome_db.chrmStartsBinCont[chr1]+start1,current_trans2,marker="o",color="b",markersize=2)
			plt.gca().set_xticks(genome_db.chrmStartsBinCont)
			plt.gca().set_xticklabels([genome_db.idx2label[i] for i in range(genome_db.chrmCount)],rotation=90,fontsize="xx-small")
#			plt.tight_layout()
			fig2.savefig(figure_path+"."+genome_db.idx2label[chr1]+"_"+str(start1*resolution/1000)+genome_db.idx2label[chr2]+"_"+str(start2*resolution/1000)+".png",dpi=600)
			plt.clf()
			plt.figure(1)


			#All trans for region1 on chr2
			current_trans = q[genome_db.chrmStartsBinCont[chr1]+start1:genome_db.chrmStartsBinCont[chr1]+end1+1,
							genome_db.chrmStartsBinCont[chr2]:genome_db.chrmEndsBinCont[chr2]].flatten()
			plt.boxplot(current_trans,positions=[len(x_labels)],whis=[1,99],showfliers=False)
			x_labels.append(genome_db.idx2label[chr1]+"_"+str(start1*resolution/1000)+"k_on "+genome_db.idx2label[chr2])

			#All trans for region2 on chr1
			current_trans = q[genome_db.chrmStartsBinCont[chr2]+start2:genome_db.chrmStartsBinCont[chr2]+end2+1,
							genome_db.chrmStartsBinCont[chr1]:genome_db.chrmEndsBinCont[chr1]].flatten()
			plt.boxplot(current_trans,positions=[len(x_labels)],whis=[1,99],showfliers=False)
			x_labels.append(genome_db.idx2label[chr2]+"_"+str(start2*resolution/1000)+"k_on "+genome_db.idx2label[chr1])
			
			#Add some interval between graph
			x_labels.append("")

	#All trans for chr1 vs chr2
	plt.boxplot(
					q[genome_db.chrmStartsBinCont[chr1]:genome_db.chrmEndsBinCont[chr1],
					genome_db.chrmStartsBinCont[chr2]:genome_db.chrmEndsBinCont[chr2]
					].flatten(),
					positions=[len(x_labels)],
					whis=[1,99],showfliers=False)
	x_labels.append(genome_db.idx2label[chr1]+"/"+genome_db.idx2label[chr2])
	x_labels.append("")
	
ax = plt.gca()
ax.set_xlim((0,len(x_labels)+1))
ax.set_xticks(range(len(x_labels)+2))
ax.set_xticklabels(x_labels, rotation=90, fontsize='xx-small')



#plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2,fontsize='xx-small')
fig1.savefig(figure_path,bbox_inches='tight',dpi=300)

#Plot chromosomal tans heatmaps for each pair of regions
# figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"Trans_"

# for i in to_plot.keys():
	# print i
	# plt.clf()
	# plotting.plot_matrix(to_plot[i])
	# plt.savefig(figure_path+i+".png",dpi=300)
	
	#Plot tans heatmaps for each pair of regions
# figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"Trans_"

#Plot tans heatmaps for each region

figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"Trans_"
# plotting.plot_matrix(t)

# ax = plt.gca()
# ax.set_xlim(0,len(q))
# xticks = [0]
# xticklables = [genome_db.idx2label[0]]
# for i in range(1,len(genome_db.chrmStartsBinCont)):
	# xticks+=range(genome_db.chrmStartsBinCont[i-1]+10,genome_db.chrmEndsBinCont[i-1],10)
	# xticklables+=(map(str,[(t-genome_db.chrmStartsBinCont[i-1])*resolution/1000 
								# for t in range(genome_db.chrmStartsBinCont[i-1]+10,genome_db.chrmEndsBinCont[i-1],10)]))
	# xticks+=[genome_db.chrmStartsBinCont[i]]
	# xticklables+=[genome_db.idx2label[i]]
	
# ax.set_xticks(xticks)
# ax.set_xticklabels(xticklables,rotation=90,fontsize="xx-small")
