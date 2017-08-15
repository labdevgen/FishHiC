import sys
try: 
	import numpy as np 
except:
	sys.path = ["/usr/lib64/python2.7/site-packages"]+sys.path
	import numpy as np
print "Numpy inported!"
from hiclib.fragmentHiC import HiCdataset
from mirnylib.systemutils import fmap,setExceptionHook
from mirnylib.genome import Genome

import os
#from defineGenome import getGenome

genomeName = "mm10"
genome_db = Genome("/mnt/storage/home/vsfishman/HiC/fasta/"+genomeName,readChrms=["#","X","Y"])

data_folder = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/"
fdataset_fname = "mESC-all-HindIII_refined.frag"
setExceptionHook()

print "Loading HiCdataset"
TR = HiCdataset(data_folder + fdataset_fname,enzymeName = "HindIII",
                        mode='r', genome=genome_db)
#print "Saving heatmap"
#TR.saveHeatmap(data_folder + fdataset_fname + "1000k.hm", 1000000)
#print "Saving heatmap"
#TR.saveHeatmap(data_folder + fdataset_fname + "40k.hm", 40000)
print "Saving heatmap"
TR.saveHeatmap(data_folder + fdataset_fname + "15k.hm", 15000,useFragmentOverlap=True)
print "Saving Hi-Res heatmap"
TR.saveHiResHeatmapWithOverlaps(data_folder + fdataset_fname + "10k_hiRes.hm", 10000)
