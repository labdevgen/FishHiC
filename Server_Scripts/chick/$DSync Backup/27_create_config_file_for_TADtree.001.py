import os
import sys
import gzip

if len(sys.argv) < 3:
	print "Usage: ",sys.argv[0]," path_to_folder_with_chromosomes_matrixes out_folder"
	print "Script will create file path_to_folder_with_chromosomes_matrixes.config in the out_folder"
	sys.exit()

in_folder_name = sys.argv[1]
out_folder = sys.argv[2]
if in_folder_name.endswith("/"):
	out_file = out_folder+"/"+in_folder_name.split("/")[-2]+".cfg"
else:
	out_file = out_folder+"/"+in_folder_name.split("/")[-1]+".cfg"
print "Writing file",out_file

files = [i for i in os.listdir(in_folder_name) if (i[0:4]=="Nchr" and i.endswith(".gz"))]
chrms=[os.path.abspath(os.path.join(in_folder_name,i)) for i in files]
lengths = [len(gzip.open(f, 'r').readline().strip().split()) for f in chrms]
labels = [i.split(".gz")[0][1:] for i in files]

chrms = [chrm for chrm,length in zip(chrms,lengths) if length > 4]
labels = [label for label,length in zip(labels,lengths) if length > 4]
lengths = [l for l in lengths if l>4]


out = open(out_file,"w")
out.write(
"""S = 50										# max. size of TAD (in bins)
M = 25										# max. number of TADs in each tad-tree
p = 3										# boundary index parameter
q = 12										# boundary index parameter
gamma = 500									# balance between boundary index and squared error in score function

""")

out.write("contact_map_path = "+",".join(chrms)+"\n")
out.write("contact_map_name = "+",".join(labels)+"\n")
out.write("N = "+",".join(map(str,[max(1,l/4) for l in lengths]))+"\n")

if not os.path.isdir(out_file+"_out"):
	os.mkdir(out_file+"_out")

out.write("output_directory = "+os.path.abspath(out_file)+"_out")