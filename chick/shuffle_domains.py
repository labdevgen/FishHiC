import numpy as np
import os
from mirnylib.systemutils import setExceptionHook
setExceptionHook()

import sys
domains_file = sys.argv[1]
data = np.genfromtxt(domains_file,dtype=np.dtype([('chrm','S10'),('start',np.int64),('end',np.int64)]),usecols = (0,1,2))
chrms = np.unique(data['chrm'])
dataByChrm={}
dataByChrmMutLengths={}
dataByChrmMutBreaks={}
dataByChrmShuffledDomains={}

for c in sorted(chrms):
	dataByChrm[c] = data[data['chrm']==c]
	lengths = dataByChrm[c]["end"] - dataByChrm[c]["start"]
	assert np.all(lengths>0)
	breaks = np.concatenate((np.atleast_1d(dataByChrm[c][0]["start"]),dataByChrm[c]["start"][1:] - dataByChrm[c]["end"][:-1]))
	assert np.all(breaks>=0)
	np.random.shuffle(lengths)
	dataByChrmMutLengths[c] = lengths
	np.random.shuffle(breaks)
	dataByChrmMutBreaks[c] = breaks
	dataByChrmShuffledDomains[c] = np.empty(shape=(len(dataByChrmMutLengths[c]),),dtype=np.dtype([('chrm','S10'),('start',np.uint32),('end',np.uint32)]))
	dataByChrmShuffledDomains[c]["chrm"] = [c]*len(dataByChrmMutLengths[c])
	dataByChrmShuffledDomains[c]["start"] = np.cumsum(dataByChrmMutLengths[c]+dataByChrmMutBreaks[c]) - dataByChrmMutLengths[c]
	dataByChrmShuffledDomains[c]["end"] = np.cumsum(dataByChrmMutLengths[c]+dataByChrmMutBreaks[c])
	assert dataByChrmShuffledDomains[c][-1]["end"] == dataByChrm[c][-1]["end"]

#print some stats on domains_file
all_lengths = data["end"]-data["start"]
print "Domains average, mean, median, std:"
print np.average(all_lengths),np.mean(all_lengths),np.median(all_lengths),np.std(all_lengths)

#save shuffled domains
for i in [""]+range(1,100):
	if os.path.isfile(domains_file+".shuffled"+str(i)):
		continue
	np.savetxt(domains_file+".shuffled"+str(i),np.hstack(tuple(dataByChrmShuffledDomains.values())),fmt="%s")
	break