f100 = "~/HiC/data/Eigensvect/ESC_full_100.eig"
f1000 = "~/HiC/data/Eigensvect/ESC_full_1000.eig"

import numpy as np

big=np.fromfile(f1000)
small=np.fromfile(f100)

res = np.array(len(big)**2)

for i in big:
	