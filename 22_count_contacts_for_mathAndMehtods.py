print "Starting programm"
from mirnylib import genome
from mirnylib import h5dict
from hiclib import binnedData
from hiclib import fragmentHiC

import numpy as np

genome_name='mm9'
genome_db = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

domain_res=1000000

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = 'Sp_full2'

heatmap_filepath=base_folder+'heatmap-res-'+str(domain_res/1000)+'KB_'+base_filename+'.hdf5'

print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(domain_res, genome_db)
BD.simpleLoad(heatmap_filepath, 'HindIII_GM_1')

q=BD.dataDict['HindIII_GM_1']

chrms = range(genome_db.chrmCount) #numner of chrms in genome
st = genome_db.chrmStartsBinCont #array of numbers of chrms start (in bins)
end = genome_db.chrmEndsBinCont #array of numbers of chrms ends (in bins). Numer end[-1] is not in chromosome (this number is 1st bin of next chrm)

def get_s(q2,N):
	s=[]
	for i in xrange(20):
		q=q2[st[i]:end[i],st[i]:end[i]]
		for ind,val in enumerate(q):
			if np.sum(val) != 0:
				s.append(np.sum(val[max(0,ind-N):min(len(val),ind+N)])/np.sum(val))
	return s

for i in range(50,100,10):
	s = get_s(q,i)
	print i,np.sum(s)/len(s)