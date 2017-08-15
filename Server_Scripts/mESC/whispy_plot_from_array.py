import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mirnylib.systemutils import setExceptionHook
setExceptionHook()
import os
dir = "/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/mm10/mESC-all-HindIII_refined.frag_res10k_hiRes.hm.IC.cis.hdf5.stat/"
files = os.listdir(dir)
files=[f for f in files if f.endswith(".dist.txt")]
arrays=dict([(file,np.loadtxt(dir+file)) for file in files])

percentiles = [1,99]
for a in arrays:
  print "min\tmax\tmean\tmedian\n"
  print np.min(arrays[a]),np.max(arrays[a]),np.mean(arrays[a]),np.median(arrays[a])
  fig, ax = plt.subplots()
  ax.set_yscale('log', basey=10)
  plt.boxplot(arrays[a],whis=percentiles,showfliers=False)
  plt.title(a+" percentiles "+"-".join(map(str,percentiles)))
  plt.savefig(dir+a+".png",dpi=300)
  plt.clf()