import numpy as np
import matplotlib.pyplot as plt

print "Creating array"
a=np.arange(10000000000).reshape(100000,100000)
print "Saving pict"
plt.matshow(a)
plt.savefig("test_mtashow.png",dpi=300)