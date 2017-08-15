print "Starting programm"
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
import numpy as np

genome_name='mm9'

genome_db1 = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
genome_db2 = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename1 = 'Sp_full'
base_filename2 = 'Fib_full'

#IMPORTANT: use iter-corrected heatmaps here. Otherwise, take care about adjustment of total reads number when calculating mask_hugeDifference
heatmap_filepath1=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename1+'.hdf5'
heatmap_filepath2=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'.hdf5'

out_heatmap_filepath2 = base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'_compressed_as_'+base_filename1+'.hdf5'
figure_path = base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'_compressed_as_'+base_filename1+'.png'

print "Loading file "+heatmap_filepath1
BD1 = binnedData.binnedData(domain_res, genome_db1)
BD1.simpleLoad(heatmap_filepath1, 'HindIII_GM_1')

print "Loading file "+heatmap_filepath2
BD2 = binnedData.binnedData(domain_res, genome_db2)
BD2.simpleLoad(heatmap_filepath2, 'HindIII_GM_1')

q1=BD1.dataDict['HindIII_GM_1']
q2=BD2.dataDict['HindIII_GM_1']


#-----DEBUG------
#print "Plotting contact matrix"
#plotting.plot_matrix(np.log(q2))
#plt.subplots_adjust(bottom=0.15)
#print "Saving figure "+figure_path+'tmp.png'
#f = open(figure_path+'tmp.png', "wb")
#plt.savefig(figure_path+'tmp.png',dpi=600)
#plt.clf()
#f.close()
#-----DEBUG END------

#DEBUG
q1=np.arange(100).reshape(10)
q2=q1*2
print q1
print q2

chrms=[0,1]
st=[0,5]
end=[5,10]
#END DEBUG

np.fill_diagonal(q1,0)
np.fill_diagonal(q2,0)


#Stage 1. Adjust total number of contacts/matrix
print "Counting total number of contacts/matrix"

#chrms = range(genome_db1.chrmCount) #numner of chrms in genome
#st = genome_db1.chrmStartsBinCont #array of numbers of chrms start (in bins)
#end = genome_db1.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)



s1 = np.sum([np.sum(q1[st[i]:end[i],st[i]:end[i]]) for i in chrms]) #total sum of contacts in 1st file (intrachr)
s2 = np.sum([np.sum(q2[st[i]:end[i],st[i]:end[i]]) for i in chrms]) #total sum of contacts in 1st file	(intrachr)

koef = float(s1)/float(s2)

#Stage 2. Adjust calculate P(s)
print "Calculating Number_of_contacts(Distance) destribution"

l=len(q1)

N1 = np.zeros(l)
N2 = np.zeros(l)

for i in chrms:
	l = genome_db1.chrmLensBin[i]
	for j in xrange(0,l):
		N1[j] = float(sum(np.diag(q1[st[i]:end[i],st[i]:end[i]],j))*2)
		N2[j] = float(sum(np.diag(q2[st[i]:end[i],st[i]:end[i]],j))*2)

N1 = N1 / s1
N2 = N2 / s2

print "Adjusting contacts according to distribution"


for i in chrms:
	l = genome_db1.chrmLensBin[i]
	for j in xrange(1,l): #Adjusting everything except main diagonal
		idx1 = range(st[i],end[i]-j) #i-indices for diagonal elements
		idx2 = range(st[i]+j,end[i]) #j-indices for diagonal elements
		if (N1[j] == 0) or (N2[j]==0):
			q2[idx1,idx2] = 0
		else:
			q2[idx1,idx2] = q2[idx1,idx2]*(N1[j]/N2[j])

for i in xrange(0,l):
	for j in xrange(i+1,l):
		q2[j,i] = q2[i,j]

print "Adjusting total contact number after correction"
s3 = float(np.sum(BD2.dataDict['HindIII_GM_1'])) / np.sum(q2)
print s3
q2=q2*s3

print q2
raise "STOP"

print "Exporting heatmap to ",out_heatmap_filepath2
BD2.dataDict['HindIII_GM_1'] = q2
BD2.export('HindIII_GM_1',out_heatmap_filepath2)

print "Plotting contact matrix"
plotting.plot_matrix(np.log(q2))
plt.subplots_adjust(bottom=0.15)
print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()