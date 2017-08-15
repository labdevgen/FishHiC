import os
genomeName = "galGal4_1MB"
enzymeName="HindIII"

#data_folder="/mnt/storage/home/vsfishman/HiC/data/chick/mapped-galGal4_1MB/"
data_folder="/mnt/storage/home/vsfishman/HiC/data/chick/mapped-galGal4_1MB/newData/"

f_out=open(data_folder+"/datasets.tsv",'w')

samples={"B":"Blood","b":"Blood","E":"ChEF","E":"ChEF"}
for i in sorted(os.listdir(data_folder)):
#  sample_name=samples[i.split("_")[0][0]]
  print i
  if not os.path.isdir(data_folder+i):
    continue
  sample_name=samples[i.split("_")[1][0]]
#  replicate="R"+i.split("_")[0][1]
  replicate="R"+i.split("_")[1][1]  
  
  for j in os.listdir(data_folder+i+"/"):
    if j.split(".")[-1] == "hdf5":
	  f_out.write(data_folder+i+"/"+j+"\t"+sample_name+"\t"+replicate+"\t"+genomeName+"\t"+enzymeName+"\n")

f_out.close()