print "Starting programm"
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Prepare to plot figure
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')


bash_output_files_dir = "/mnt/storage/home/vsfishman/HiC/data/BASH/bash_output/Fib_full2_100KB_all_domains.txt/"

f_res = open(bash_output_files_dir+"/hdratios.txt","w")
f_res.write("Sample\tChr\tStart\tEnd\tHD-ration\tNumber_of_Bins\tNormalized_HDration\n")

d = [i for i in os.listdir(bash_output_files_dir) if os.path.isdir(bash_output_files_dir +"/"+i)]
count = 0

for i in sorted(d):
	if (count % 100 == 0):
		print "Processed",count
	count += 1
	dir = bash_output_files_dir + "/" + i + "/"
	if not os.path.isfile(dir + "mode_p.txt"):
		print "WARNING! FILE ",dir + "mode_p.txt"," NOT FOUND"
		continue
	with open(dir + "mode_p.txt") as f_in:
		m=[]
		for line in f_in:
			x=float(line.split("\t")[0])
			y=float(line.split("\t")[1])
			z=float(line.split("\t")[2])
			m.append([x,y,z])
	m=np.array(m)
	for q in xrange(0,3): #changing coordinate x,y,z, to x-<x>,y-<y>,z-z<z>. Bringing wight center to 0,0,0 position
		m[:,q] -= np.average(m[:,q])
	#Plot before PCA
#	ax.scatter(m[:,0], m[:,1], m[:,2],c='r',marker='^') #plotting figure 3D
	Y=m/(len(m[:,0])**0.5)
	u,s,v=np.linalg.svd(Y) #Singular value decomposition of m
	m_new = np.dot(m,v) #getting new coordinates
	
	#Plot after PCA
#	ax.scatter(m_new[:,0], m_new[:,1], m_new[:,2],c='b',marker='o') #plotting figure 3D
#	ax.scatter(0, 0, 0,c='y',marker='*') #plotting dot in the 0,0,0 point
	
	#Set axies labels
#	ax.set_xlabel('x') 
#	ax.set_ylabel('y')
#	ax.set_zlabel('z')

	height = abs(np.percentile(m_new[:,0],90)-np.percentile(m_new[:,0],10))
	diameter = [(m_new[q,1]**2+m_new[q,2]**2)**0.5 for q in xrange(len(m_new[:,1]))]
	diameter = np.average(diameter)
	hd_ratio = height/diameter
	
#	f = open(dir+"picure.png", "wb")
#	plt.savefig(dir+"picure.png",dpi=300)
#	f.close()
	temp = i.split(".")
	f_res.write(temp[-4] + "\t" + temp[-3] + "\t" +temp[-2] + "\t" +temp[-1] + "\t" + str(hd_ratio) + "\t" + str(len(m_new[:,0])))
	hd_ratio_normalized = hd_ratio/len(m_new[:,0]) #hd-ration normizlized to the number of genomic regions
	f_res.write("\t"+str(hd_ratio_normalized) + "\n")
print "Done"