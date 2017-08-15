import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lib_matrix_operations import *
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

from mirnylib.numutils import PCA, EIG, correct, \
    ultracorrectSymmetricWithVector, isInteger, \
    observedOverExpected, ultracorrect, adaptiveSmoothing, \
    removeDiagonals

genome_name='mm9'
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X','Y'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'Sp_full2_Y'

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5-raw'
#figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb-intrarchr-new4.png'
figure_path=base_folder+base_filename+"_"+str(domain_res/1000)+'kb-intrarchr-newN.png'

genome_fai_filepath='../fasta/'+genome_name+'/'+genome_name+'.fai'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')
q=BD.dataDict['HindIII_GM_1']
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

		print i,j,p1*(total/2.0),p2*(total/2.0),q2[i][j], znam
		q2[i][j] = q2[i][j]/znam
print q2

q1=np.array(q2)
for i in xrange(len(q1)):
	q1[i,i] = 1.0
	
#q1=np.log2(q1)
print q1

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
plotting.plot_matrix(q1,cmap='seismic')
plt.subplots_adjust(bottom=0.15)
print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()
