print "Starting programm"
import os

bash_input_files_dir = "/mnt/storage/home/vsfishman/HiC/data/BASH/bash_input/Sp_full2_100KB_all_domains.txt/"
q_files_dir = "/mnt/storage/home/vsfishman/HiC/data/BASH/bash_q/"
bash_output_files_dir = "/mnt/storage/home/vsfishman/HiC/data/BASH/bash_output/"

postfix = bash_input_files_dir.split("/")
if postfix[-1] != "":
	postfix=postfix[-1]
else:
	postfix=postfix[-2]

bash_output_files_dir += "/"+postfix
q_files_dir_full = q_files_dir + "/"+postfix
if not os.path.isdir(q_files_dir_full):
	os.mkdir(q_files_dir_full)



if not os.path.isdir(bash_output_files_dir):
	os.mkdir(bash_output_files_dir)

f=os.listdir(bash_input_files_dir)

number_of_files = len(f)
jobs_pro_file = max(number_of_files/135,1)

def print_header(f_o):
	f_o.write("#!/bin/bash \n#PBS -V -r n \n#PBS -l select=1:ncpus=6:mem=12gb,cput=64:10:00,walltime=32:15:00\n#PBS -q bl2x220g7q\n#PBS -N BASH\n#PBS -j oe\ndate\n#!/bin/sh\ndate\n")
	f_o.write("cd $HOME/HiC/bin/BACH_src/\n")

open_new_f=True
j = 0
count = 0
for i in sorted(f):
	if i.split(".")[-1] == "heatmap":
		count += 1
		h = bash_input_files_dir+"/"+i
		base = ""
		for q in i.split(".")[:-1]:
			base += q+"."
		base = base[:-1]
		c = bash_input_files_dir+"/"+base+".feature"
		
		if open_new_f:
			f_out = open(q_files_dir_full + "/"+i+"_and_"+str(jobs_pro_file)+"next_jobs.q","w")
			print_header(f_out)
			open_new_f = False
		f_out.write("if [ -e "+bash_output_files_dir+"/"+base+"/mode_p.txt ];\nthen echo \"file exists...skipping!\";\nelse ")
		f_out.write("mkdir -p "+bash_output_files_dir+"/"+base+"\n")
		f_out.write("./BACH  -i "+h+" -v "+c +" -K 100 -MP 10 -NG 5000 -NT 50 -L 50 -SEED 1 -o "+bash_output_files_dir+"/"+base+"\n")
		f_out.write("fi;")
		f_out.write("echo \"Run complted\"\n")
		j += 1
		
		if j>=jobs_pro_file:
			f_out.close()
			open_new_f = True
			j=0
print count," preoceeded"
print "Done!"
