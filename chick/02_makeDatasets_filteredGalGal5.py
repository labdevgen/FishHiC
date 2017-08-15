import os
genomeName = "GalGal5filtered"
enzymeName="HindIII"

#data_folder="/mnt/storage/home/vsfishman/HiC/data/chick/mapped-galGal4_1MB/"
data_folder="/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/chick/mapped-GalGal5filtered/"

f_out=open(data_folder+"/datasets.tsv",'w')

samples={"B":"Blood","E":"ChEF"}
for i in sorted(os.listdir(data_folder)):
#  sample_name=samples[i.split("_")[0][0]]
  print i
  if not os.path.isdir(data_folder+i):
    continue
  if i[0:4] == "exp2":
    sample_name=samples[i.split("_")[2][0]]
    replicate = i.split("_")[2][1]
  else:
    sample_name=samples[i.split("_")[1][0]]
    replicate = i.split("_")[1][1]

	
  replicate="R"+replicate
  
  for j in os.listdir(data_folder+i+"/"):
    if j.split(".")[-1] == "hdf5":
	  f_out.write(data_folder+i+"/"+j+"\t"+sample_name+"\t"+replicate+"\t"+genomeName+"\t"+enzymeName+"\n")

f_out.close()