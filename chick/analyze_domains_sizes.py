import numpy as np
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab


domains = {
	"Dixon_Blood_40kb":"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsBlood_all_HindIII_40k.hm.IC_domains_40KB/",
	"Dixon_ChEF_40kb":"/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomainsChEF_all_HindIII_40k.hm.IC_domains_40KB/",
	"Armatus_Blood_25kb":"mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-25k.hm.gzipped_matrix/",
	"Armatus_Blood_40kb":"mapped-GalGal5filtered/GalGal5filteredChrmLevel/Blood-all-HindIII-40k.hm.gzipped_matrix/",
	"Armatus_ChEF_25kb":"mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-25k.hm.gzipped_matrix/",
	"Armatus_ChEF_40kb":"mapped-GalGal5filtered/GalGal5filteredChrmLevel/ChEF-all-HindIII-40k.hm.gzipped_matrix/"
}

percentiles = [1,99]
fig, ax = plt.subplots()
ax.set_yscale('log', basey=10)

def add_graph(files):
	sizes = []
	for f in files:
		if len(open(f).readlines()) == 0:
			continue
		data = np.genfromtxt(f,dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]),usecols = (0,1,2))
		tmp_sizes = data["end"]-data["start"]
		assert np.all(tmp_sizes>0)
		sizes.append(np.atleast_1d(tmp_sizes))
	sizes = np.concatenate(tuple(sizes))
	print len(sizes),sum(sizes),sum(sizes)*100/1021439028.,np.average(sizes),np.median(sizes),sizes.min(),sizes.max(),np.std(sizes),np.mean(sizes), \
			np.percentile(sizes,percentiles[0]),np.percentile(sizes,percentiles[1]),np.percentile(sizes,25),np.percentile(sizes,75)
	return sizes

counter = 1
plots = []
cmap_blue = matplotlib.cm.get_cmap('Blues')
cmap_green = matplotlib.cm.get_cmap('Greens')

for d in sorted(domains):
	print d
	
	if "Dixon" in d:
		files=[domains[d]+f for f in os.listdir(domains[d]) if f.endswith(".final_domains")]
		data = add_graph(files)
		bp = plt.boxplot([data],positions=[counter],whis=percentiles,showfliers=False,patch_artist=True)
		if "Blood" in d:
			pylab.setp(bp['boxes'], color='blue')
		else:
			pylab.setp(bp['boxes'], color='green')
		pylab.setp(bp['boxes'],hatch = '+')
		counter += 1
	if "Armatus" in d:
		files=[f for f in os.listdir(domains[d]) if ".domains." in f]
		gamma_options = [f.split(".domains.")[-1].split(".txt")[0] for f in files]
		gamma_options = list(set(gamma_options))
		for ind_g,g in enumerate(sorted(gamma_options)[::1]):
			print g
			data=add_graph([domains[d]+f for f in files if g in f])
			bp = plt.boxplot([data],positions=[counter],whis=percentiles,showfliers=False,patch_artist=True)
			if "Blood" in d:
				pylab.setp(bp['boxes'], color=cmap_blue(float(ind_g)/len(gamma_options[::1])))
			else:
				pylab.setp(bp['boxes'], color=cmap_green(float(ind_g)/len(gamma_options[::1])))
			counter += 1
ax.set_xlim([0,counter+1])

plt.savefig("domains_stats.png",dpi=300)
