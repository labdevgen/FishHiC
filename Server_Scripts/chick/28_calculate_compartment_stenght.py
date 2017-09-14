#Please note
#This code is mainly based on the doSaddle function from singleShared.py
#Which is a part of Mirny Lab HiCLib libraries, version 85979ac
import sys
import numpy as np
import math
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from mirnylib.numutils import observedOverExpected
from mirnylib.systemutils import setExceptionHook
setExceptionHook()



if len(sys.argv) < 3:
	print "Usage: ",sys.argv[0]+" heatmap_file E1_file"
	sys.exit()

heatmap_filepath = sys.argv[1]
E1_file = sys.argv[2]

sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
import figPath
import ntpath
figure_path=figPath.figure_path+ntpath.basename(heatmap_filepath)+"_"+'Compartment_strength'
if not os.path.isdir(figure_path):
	os.mkdir(figure_path)

figure_path+="/"

print "Saving results to",figure_path


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

def doSaddles(q,E1_values,genome_db):
	saddles={}
	for chrom in range(genome_db.chrmCount):
		saddle = np.ones((5,5), dtype = float)
		st = genome_db.chrmStartsBinCont[chrom]
		end = genome_db.chrmEndsBinCont[chrom]
		cur = q[st:end,st:end]
		
		E1 = E1_values[st:end]

		mask = np.sum(cur , axis=0) > 0
		if sum(mask)>5:
				cur = cur [mask]
				cur = cur [:,mask]
				
				cur = observedOverExpected(cur)
				E1 = E1[mask]
				assert cur.shape[0]==cur.shape[1]==len(E1)
				
				for i in range(5):
					for j in range(5):
						P1, P2 = np.percentile(E1, [20 * i, 20 * i + 20])
						mask1 = (E1 > P1) * (E1 < P2)
						P1, P2 = np.percentile(E1, [20 * j, 20 * j + 20])
						mask2 = (E1 > P1) * (E1 < P2)
						if sum(mask1)*sum(mask2) != 0:
							saddle[i, j] = np.nanmean(cur[np.ix_(mask1, mask2)])
						else:
							saddle[i, j] = None
				saddles[genome_db.idx2label[chrom]]=saddle
		else:
			pass
			#print "Ommiting chromsome ",genome_db.idx2label[chrom]
	
	all_average=np.zeros((5,5),dtype=float)
	for i in range(5):
		for j in range(5):
			all_average[i,j] = np.average([saddles[c][i,j] for c in saddles if not np.isnan(saddles[c][i,j])])
	saddles["all_average"] = all_average
	strength = math.log(all_average[0,0]*all_average[-1,-1]/(all_average[0,-1]*all_average[0,-1]))
	return saddles,strength

#Load hmap
resolution = int(heatmap_filepath.split("-")[-1].split("k")[0])*1000
print "Resolution determined: ",resolution
print "Loading file "+heatmap_filepath
BD = binnedData.binnedData(resolution, genome_db)
BD.simpleLoad(heatmap_filepath, 'heatmap')
q=BD.dataDict['heatmap']

#Load E1
E1_values = np.genfromtxt(E1_file,dtype=None)['f2']
assert len(E1_values)==len(q)

saddles,strength = doSaddles(q,E1_values,genome_db)

vmin = np.min(saddles["all_average"])
vmax = np.max(saddles["all_average"])
print "Min = ",vmin,"Max = ",vmax

#Shufle E1 to get boostrap control
import pickle
SD=[]
for i in range(100):
	if i%5 == 0:
		print "Bootstraping ",i,"%"
	np.random.shuffle(E1_values)
	SD.append(doSaddles(q,E1_values,genome_db))
with open(figure_path+"bootstrup_dump","w") as f:
	pickle.dump([SD],f)

print "Strength=",strength,"+/-",np.std([i[1] for i in SD])
print "Bootstrap average = ",np.average([i[1] for i in SD])


for i in saddles:
	print "saving ",i
	plt.clf()
	np.savetxt(figure_path+i+".txt",saddles[i])
	plotting.plot_matrix(np.log(saddles[i]))
	if i=="all_average":
		plt.title("Strength="+str(strength))
	plt.savefig(figure_path+i+".png",dpi=300)