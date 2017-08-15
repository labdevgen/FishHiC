print "Starting programm"
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
length_step = 20000000 # length is in KB

if ((length_step % domain_res != 0) or (domain_res >= length_step)):
	print "Length step ", length_step," incompatible with resolution ",domain_res
	print "length_step % domain_res =", (length_step % domain_res)
	raise "ERROR!"
else:
	length_step_inBin = length_step/domain_res
	print "Length of distance step is ",length_step_inBin," bins"

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename1 = 'Sp_full'
base_filename2 = 'Fib_full'

#IMPORTANT: use iter-corrected heatmaps here. Otherwise, take care about adjustment of total reads number when calculating mask_hugeDifference
heatmap_filepath1=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename1+'.hdf5'
heatmap_filepath2=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename2+'.hdf5'

use_weights = 2 #mask all points where difference is significant, but <= use_weights times. 1 - for using no weights
if (use_weights == 1):
	figure_path=base_folder+base_filename1+"_vs_"+base_filename2+"_"+str(domain_res/1000)+"kb-stats_diff_Step"+str(length_step_inBin)+"bins.png"
else:
	figure_path=base_folder+base_filename1+"_vs_"+base_filename2+"_"+str(domain_res/1000)+"kb-stats_diff_Step"+str(length_step_inBin)+"bins"+str(use_weights)+"timesDif.png"
	print "WARNING: Assuming datasets iteratively corrected"

print "Loading file "+heatmap_filepath1
BD1 = binnedData.binnedData(domain_res, genome_db1)
BD1.simpleLoad(heatmap_filepath1, 'HindIII_GM_1')

print "Loading file "+heatmap_filepath2
BD2 = binnedData.binnedData(domain_res, genome_db2)
BD2.simpleLoad(heatmap_filepath2, 'HindIII_GM_1')

q1=BD1.dataDict['HindIII_GM_1']
q2=BD2.dataDict['HindIII_GM_1']

if (use_weights != 1):
	print "Adjusting total N of contacts/bin"
	for i in q1:
		s1 = sum(i)
		if s1 != 0:
			break

	for i in q2:
		s2 = sum(i)
		if s2 != 0:
			break

#q2=np.array(q1)
#from random import randint
#for i in xrange(100):
#	for j in xrange(i,100):
#		q2[i,j]+=randint(0,100)
#		q2[j,i] = q2[i,j]

print "Calcultaing total number of interactions"

distr=scipy.stats.norm(0,1) #normal distribution with disp,E

print "Setting mapability mask"
mapability = 2 #mapability is a minmal number of reads connecting 2 points to consider these points in future
if mapability != 1:
	print "WARNING: mapability calculated from IC-heatmaps. Absolute values of contact numbers used, however, relative numbers are used in IC-heatmaps"
#masks are arrays where in positions that I need to remove are 1(True)
mask_mapability = np.logical_or((q1 < mapability), (q2 < mapability)) #all contacts with >=0 reads

if use_weights != 1:
	print "Setting weights mask"
	koeff = float(s2)/s1
	mask_hugeDifference = abs(np.log10((q1*koeff)/q2)) < math.log10(use_weights) 

print "Calculating Number_of_contacts(Distance) destribution"

N1 = np.ones_like(q1) #in this array total number of contact at distance x will be stored
N2 = np.ones_like(q2) #in this array total number of contact at distance x will be stored
l=len(N1)

np.fill_diagonal(q1,0)
np.fill_diagonal(q2,0)

for i in xrange(0,l,length_step_inBin):
	idx1=[] #array of diagonal indexes to modify
	idx2=[] #array of diagonal indexes to modify
	for j in xrange(i,i+length_step_inBin):
		idx1 += range(l-j) #i-indices for diagonal elements
		idx2 += range(j,l) #j-indices for diagonal elements
	N1[idx1,idx2] = sum([sum(np.diag(q1,k)) for k in xrange(i,i+length_step_inBin)])
	N2[idx1,idx2] = sum([sum(np.diag(q2,k)) for k in xrange(i,i+length_step_inBin)])

if ((l % length_step_inBin) != 0):
	last_it = (l/length_step_inBin)*length_step_inBin
	idx1=[] #array of diagonal indexes to modify
	idx2=[] #array of diagonal indexes to modify
	for j in xrange(last_it,l):
		idx1 += range(l-j)
		idx2 += range(j,l)
	N1[idx1,idx2] = sum([sum(np.diag(q1,k)) for k in range(last_it,l)])
	N1[idx1,idx2] = sum([sum(np.diag(q2,k)) for k in range(last_it,l)])
	
#__WARNING: q1 AND q2 ARE MODIFIED HERE. DO NOT USE THESE VARIABLES LATER!!!____
q1 += mask_mapability #This will remove zeros from q1 and q2 
q2 += mask_mapability #This will remove zeros from q1 and q2

print "Calculating probabilities"
p1=q1/N1
p2=q2/N2

print "Calculating Q-statistic"
Q=abs((p1-p2)/((p1*(1-p1)/N1+p2*(1-p2)/N2)**0.5))
print "Calculating p-values"

number_of_nonzero_ineractions = len(mask_mapability)**2-np.sum(mask_mapability)
count = 0
step = len(Q)/10

for i in xrange(len(Q)):
	if (i % step) == 0:
		print (i/step)*10,"%"
	for j in xrange(i,len(Q)):
		if mask_mapability[i,j]:
			Q[i,j]=1 #This points will be masked later anyway
			Q[j,i]=1 #This points will be masked later anyway			 
			continue
		Q[i,j]=2-2*distr.cdf(Q[i,j])
		Q[i,j] *= number_of_nonzero_ineractions #Calculating p-value as p[i,j]*N(nuber of examined contacts)
		if Q[i,j] > 1:
			Q[i,j]=1 # make all p-valut >1 equal to 1
		Q[j,i]=Q[i,j]
		count += 1

print "Counted ", count ," p-values"

print "number_of_nonzero_ineractions =",number_of_nonzero_ineractions," of ",len(mask_mapability)**2,"(",(float(number_of_nonzero_ineractions)/len(mask_mapability)**2)*100,"%)"

alpha=0.05 #level of significance

mask_pval = np.logical_and(np.logical_not(mask_mapability), Q>alpha) #all contacts p-value >=0.05 reads and they are not in mask

N_significant = number_of_nonzero_ineractions-np.sum(mask_pval)
print "Number of values with p-val <= ",alpha," = ",N_significant ," out of ",number_of_nonzero_ineractions, "(",(float(N_significant)/number_of_nonzero_ineractions)*100,"%)"

if use_weights != 1:
	mask_hugeDifference = np.logical_and(np.logical_not(mask_mapability), mask_hugeDifference) # take out elements that are already in mapability
	mask_hugeDifference = np.logical_and(np.logical_not(mask_pval), mask_hugeDifference) # take out elements that are already in p-val
	N_hugeDif = N_significant-np.sum(mask_hugeDifference)
	print "Number of values with difference >= ",use_weights," times = ",N_hugeDif," out of ",N_significant, "(",(float(N_hugeDif)/N_significant)*100,"%)"

print "Applying mask_mapability"
Q = Q*np.logical_not(mask_mapability)+mask_mapability*1000.0 #all points where reads were not mapable set to 1000
print "Applying p-value mask"
Q = Q*np.logical_not(mask_pval)-mask_pval*500.0 #all points where reads were p-value > alpha set to -1000
if use_weights != 1:
	print "Applying weights mask"
	Q = Q*np.logical_not(mask_hugeDifference)+mask_hugeDifference*500.0 #all points where difference is less then use_weights times

if domain_res < 1000000:
	Q = Q[0:5000,0:5000]
	print "Warning: saving only a part of picture due to high resolution"

print "Plotting contact matrix"
plotting.plot_matrix(Q)
plt.subplots_adjust(bottom=0.15)
print "Saving figure "+figure_path
f = open(figure_path, "wb")
plt.savefig(figure_path,dpi=600)
f.close()