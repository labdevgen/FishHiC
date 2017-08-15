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

genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

domain_res=1000000


heatmap_filepath=sys.argv[1]
#figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb-intrarchr-new4.png'
figure_path=heatmap_filepath+"_"+str(domain_res/1000)+'kb-intrarchr-newN.png'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'heatmap')
q=BD.dataDict['heatmap']
q1=[]
print "Calculating matrix"

#genome_db.chrmCount=3
for i in range(genome_db.chrmCount):
	print "Calculating chr ",i+1
	q1.append([])
	for j in range(genome_db.chrmCount):
		#remove diagonal elements
		if i == j:
			q1[i].append(0.0)
			continue
		chr_i_start_bin=sum(genome_db.chrmLensBin[0:i])
		chr_i_end_bin=sum(genome_db.chrmLensBin[0:i+1])

		chr_j_start_bin=sum(genome_db.chrmLensBin[0:j])
		chr_j_end_bin=sum(genome_db.chrmLensBin[0:j+1])

		s=0.0
		for k in range(chr_i_start_bin,chr_i_end_bin):
			for n in range(chr_j_start_bin,chr_j_end_bin):
				s += q[k,n]
#		s=s/(genome_db.chrmLensBin[i]*genome_db.chrmLensBin[j])
#		s=s/100000.0
		q1[i].append(s)

total=sum([i for i in [sum(j) for j in q1]])

print "chr_i\tchr_j\tchr_i_length\tchr_j_length\ttcontacts_between"

for i in range(genome_db.chrmCount):
	for j in range(i+1,genome_db.chrmCount):
		print i,j,genome_db.chrmLensBin[i],genome_db.chrmLensBin[j],q1[i][j]

for i in q1:
	s=""
	for j in i:
		s+=str(j)+"\t"
	print s, sum(i)

q1=np.array(q1)
q2=np.array(q1)



#print q2
for i in range(genome_db.chrmCount):
	for j in range(genome_db.chrmCount):
		if i==j: continue
		s1=float(sum(q1[i]))
		s2=float(sum(q1[j]))
		p1=(s1/total)*(s2/(total-s1))
		p2=(s2/total)*(s1/(total-s2))

#		p1=(s1/total)*(s2/(total))
#		p2=(s2/total)*(s1/(total))

		znam = (p1+p2)*(total/2.0)

		print genome_db.idx2label[i],genome_db.idx2label[j],p1*(total/2.0),p2*(total/2.0),q2[i][j], znam,s1,s2
		q2[i][j] = q2[i][j]/znam
print q2

q1=np.array(q2)
for i in xrange(len(q1)):
	q1[i,i] = 1.0
	
q1=np.log2(q1)
#print q1

#currentEIG, eigenvalues = EIG(q1, numPCs=3)
#if currentEIG[0][-1] > currentEIG[0][-2]:
#	currentEIG[0]=-currentEIG[0]
#print currentEIG

#print "Saving figure "+figure_path+".eig.png"
#f = open(figure_path+".eig.png", "wb")
#plt.plot(range(len(currentEIG[0])),currentEIG[0],"bo")
#plt.savefig(figure_path+".eig.png",dpi=600)
#f.close()
#plt.clf()

#f1 = open(figure_path+".eig", "wb")
#for i in range(len(currentEIG[0])):
#	s = "chr"+str(i)+"\t"+str(currentEIG[0][i])+"\n"
#	f1.write(s)
#f1.close()

#mi=0.8 #Min and Max values for pictures
#ma=1.2 #To make colors same
#if (np.max(q1) > ma) or (np.min(q1) < mi):
#	print "Current max and min are ",np.max(q1),np.min(q1)
#	raise Exception("Array values out of "+str(mi)+"\t"+str(ma))

print "Plotting contact matrix"
#plotting.plot_matrix(q1,vmin=mi, vmax=ma, cmap='seismic')
q3 = np.zeros_like(q2)
order = np.array([genome_db.idx2label[i] for i in range(genome_db.chrmCount)])
def sortfunct(s):
		if s in ["chrW","chrZ"]:
			return 10000000
		else:
			return 100*len(s)

sortorder = sorted(range(len(order)), key=lambda k: (sortfunct(order[k]),order[k]))
print order
print sortorder
print order[sortorder]
for i in xrange(len(q2)):
	for j in xrange(len(q2)):
		q3[i,j] = q2[sortorder[i],sortorder[j]]

print order[sortorder]
plotting.plot_matrix(q3,cmap='seismic',vmin=-1.5,vmax=4,ticklabels1=order[sortorder])
plt.subplots_adjust(bottom=0.15)
print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=300)
f.close()

for ind,st,end in zip(sortorder,genome_db.chrmStartsBinCont[sortorder],genome_db.chrmEndsBinCont[sortorder]):
	L = genome_db.chrmLens[ind]
	Total_L = sum(genome_db.chrmLens)
	intra = float(np.sum(q[st:end,st:end]))
	inter = np.sum(q[st:end])-intra
	intra /= 2.
	
	intra = intra / L**2
	inter = inter / (L*(Total_L-L))
	print genome_db.idx2label[ind],intra/inter