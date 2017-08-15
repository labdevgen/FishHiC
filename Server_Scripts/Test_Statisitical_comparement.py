import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lib_matrix_operations import *
import scipy.stats
import os

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC
import math

genome_name='mm9'

genome_db1 = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])
genome_db2 = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename1 = 'Sp_full'
base_filename2 = 'Fib_full'

heatmap_filepath1=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename1+'.hdf5'
heatmap_filepath2=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'.hdf5'

figure_path=base_folder+base_filename1+"_vs_"+base_filename2+"_"+str(domain_res/1000)+'kb-stats_diff.png'

print "Loading file "+heatmap_filepath1
BD1 = binnedData.binnedData(domain_res, genome_db1)
BD1.simpleLoad(heatmap_filepath1, 'HindIII_GM_1')

print "Loading file "+heatmap_filepath2
BD2 = binnedData.binnedData(domain_res, genome_db2)
BD2.simpleLoad(heatmap_filepath2, 'HindIII_GM_1')

q1=BD1.dataDict['HindIII_GM_1']
#q2=BD2.dataDict['HindIII_GM_1']

#calcultaing test array
print "calcultaing test array"
q1=np.zeros((6,6))
q1 += 100.0
q1[0,0] = 0
q2=np.array(q1)


from random import randint
for i in xrange(3):
	for j in xrange(i,3):
		q2[i,j] += randint(-50,50)
		q2[j,i] = q2[i,j]

print q1
print q2
#for i,iv in enumerate(q2):
#	for j,jv in enumerate(iv):
#		if ((i % 10) in range(0,5)) or ((j % 10) in range(0,5)):
#			q2[i,j] = q2[i,j]*100

print "Calcultaing 1"
n1 = np.sum(q1)/2.0 # total number of interactions in 1st file
n2 = np.sum(q2)/2.0 # total number of interactions in 2nd file

distr=scipy.stats.norm(0,1) #normal distribution with disp,E

print "Calcultaing 2"

p1=q1/n1
p2=q2/n2
Q=abs((p1-p2)/((p1*(1-p1)/n1+p2*(1-p2)/n2)**0.5))
Q=2-2*distr.cdf(Q)
print "After Q=2-2*distr.cdf(Q)"
print Q
print Q < 0.05
#mask is array where in positions that I need to remove are 1
mapability = 1
mask = np.logical_or((q1 <= mapability), (q2 <= mapability)) #all contacts with >=0 reads

number_of_nonzero_ineractions = np.sum(mask)

Q = Q*number_of_nonzero_ineractions

mask2 = np.logical_and(np.logical_not(mask), Q>=0.05) #all contacts p-value >=0.05 reads and they are not in mask

def mymin(x):
	return min(x,1)

array_elwise_min = np.vectorize(mymin)

print np.sum(mask)
print np.sum(mask2)

Q = array_elwise_min(Q)

print "mask"
print mask

print "mask2"
print mask2

Q = Q*np.logical_not(mask)+mask*1000 
Q = Q*np.logical_not(mask2)-mask2*1000

print Q

#Q = q1

#l = len(q1)
#for i in xrange(l):
#	for j in xrange(i,l):
#		if (q1[i,j]<=0 or q2[i,j]<=0):
#			Q[i,j]= -1
#			number_of_zero_ineractions += 1
#			continue
#		else:
#			p1 = q1[i,j]/n1
#			p2 = q2[i,j]/n2
#			Q[i,j] = abs((p1-p2)/((p1*(1-p1)/n1+p2*(1-p2)/n2)**0.5))
#			Q[i,j] = 2-2*distr.cdf(Q[i,j])
#
#print "Ajusting p-values"
#l2 = (l**2)/2-number_of_zero_ineractions
#for i in xrange(l):
#	for j in xrange(i,l):
#		if (Q[i,j] != -1):
#			Q[i,j] = min(1,Q[i,j]*l2)
#			if (Q[i,j] > 0.05):
#				Q[i,j] = 1000
#		Q[j,i]=Q[i,j]

#Q = Q*((len(Q)**2)-number_of_zero_ineractions)

#print "Plotting contact matrix"
#plotting.plot_matrix(Q)
#plt.subplots_adjust(bottom=0.15)
#print "Saving figure "+figure_path
#f = open(figure_path, "wb")
#plt.savefig(figure_path,dpi=600)
#f.close()