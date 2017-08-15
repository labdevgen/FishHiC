import os
genomeName = "mm10"

Hind_III_dsets = ["SRR"+str(i) for i in range(400251,400257)+range(443884,443886)+range(1504836,1504839)]
NcoI_dsets = ["SRR"+str(i) for i in range(443886,443889)+range(1504839,1504841)]

#data_folder="/mnt/storage/home/vsfishman/HiC/data/chick/mapped-galGal4_1MB/"
data_folder="/mnt/storage/home/vsfishman/HiC/data/mESC/mapped-mm10/"

f_out=open(data_folder+"/datasets.tsv",'w')

count = []
for i in sorted(os.listdir(data_folder)):
  if not os.path.isdir(data_folder+i):
    continue
#  sample_name=samples[i.split("_")[0][0]]
  sample_name="mESC"
#  replicate="R"+i.split("_")[0][1]
  replicate="R1"
  if i in Hind_III_dsets:
    enzymeName = "HindIII"
  elif i in NcoI_dsets:
    enzymeName = "NcoI"
  else:
    print "!!!!!!!",i
    raise
  count += [i]
  for j in sorted(os.listdir(data_folder+i+"/")):
    if j.split(".")[-1] == "hdf5":
	  f_out.write(data_folder+i+"/"+j+"\t"+sample_name+"\t"+replicate+"\t"+genomeName+"\t"+enzymeName+"\n")

if len(count) != len(Hind_III_dsets+NcoI_dsets):
	print count,len(Hind_III_dsets+NcoI_dsets)
	print [i for i in Hind_III_dsets+NcoI_dsets if not i in count]
	raise
f_out.close()
