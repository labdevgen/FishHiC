import numpy as np
import scipy

from mirnylib import genome
from mirnylib import h5dict
from hiclib import binnedData

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

colors=['b','g','r','c','m','y','k','w']
markers=['o','v','s','p','*','+','x','D']

inner_window_size=None

def compare_nonzerro(array1,array2):
  res=0
  for i in range(len(array1)):
	if ((array1[i]==0 and array2[i]==0) or (array1[i]<>0 and array2[i]<>0)):
		res += 1
  
  return [float(res)/float(len(array1)),0]

def window_sperman(array1,array2,window_size=2,list_analize_type="Sp"):
  cur_pos=0
  result = []
  temp = []
  if hasattr(array1[0], '__iter__'):
		for i in range(0,len(array1)):
			s1=sum(array1[i]) #Calculate sum in line
			s2=sum(array2[i]) #For each array
			if (s1==0 or s2==0): # empty line in contact matrix
				temp.append(['Undef','Undef']) #!!!!!!!!!!!!!!!!!!!-Change it later-!!!!!!!!!!!!!
				continue
			#STEP I - analize current element of array and put result into temp variable
			if len(list_analize_type.split("_w")) > 1:
				if (inner_window_size==None):
					global inner_window_size
					inner_window_size=int(raw_input("Enter inner Window size ").strip())
				temp_result=window_sperman(array1[i],array2[i],inner_window_size,list_analize_type.split("_w")[0])
				mean=np.mean([row[0] for row in temp_result])
				max_err=max([row[1] for row in temp_result])
				temp.append([mean,max_err])
			elif list_analize_type=="Sp":
				temp.append(scipy.stats.spearmanr(array1[i],array2[i]))
			elif list_analize_type=="pears":
				temp.append(scipy.stats.pearsonr(array1[i],array2[i]))
			elif list_analize_type=="Evkl":
				temp.append([np.linalg.norm(np.array([a/s1 for a in array1[i]])-np.array([b/s2 for b in array2[i]])),0])
			elif list_analize_type=="proc":
				temp.append(compare_nonzerro(array1[i],array2[i]))
			else:
				raise "Unknown list_analize_type"
		#STEP II - Split resulted temp array into windows and make an averege throgh the window
		while (cur_pos+window_size < len(temp)):
			if 'Undef' in [row[0] for row in temp[cur_pos:cur_pos+window_size]]:
				result.append(['Undef','Undef']) #!!!!!!!!!!!!!!!!!!!-Change it later-!!!!!!!!!!!!!
			else:
				mean=np.mean([row[0] for row in temp[cur_pos:cur_pos+window_size]])
				max_err=max([row[1] for row in temp[cur_pos:cur_pos+window_size]])
				result.append([mean,max_err])
			cur_pos = cur_pos+window_size
		if 'Undef' in [row[0] for row in temp[cur_pos:]]:
			result.append(['Undef','Undef']) #!!!!!!!!!!!!!!!!!!!-Change it later-!!!!!!!!!!!!!
		else:
			mean=np.mean([row[0] for row in temp[cur_pos:]])
			max_err=max([row[1] for row in temp[cur_pos:]])
			result.append([mean,max_err])
  else:
	while (cur_pos+window_size < len(array1)):
		if list_analize_type=="Sp":
			result.append(scipy.stats.spearmanr(array1[cur_pos:cur_pos+window_size-1],array2[cur_pos:cur_pos+window_size-1]))
		elif list_analize_type=="pears":
			temp.append(scipy.stats.pearsonr(array1[i],array2[i]))
		elif list_analize_type=="Evkl":
			result.append([np.linalg.norm(np.array([a/s1 for a in array1[cur_pos:cur_pos+window_size-1]])-np.array([b/s2 for b in array2[cur_pos:cur_pos+window_size-1]])),0])
		elif list_analize_type=="proc":
			result.append(compare_nonzerro(array1[cur_pos:cur_pos+window_size-1],array2[cur_pos:cur_pos+window_size-1]))
		else:
			raise "Unknown list_analize_type"
		cur_pos = cur_pos+window_size
	if list_analize_type=="Sp":
		result.append(scipy.stats.spearmanr(array1[cur_pos:],array2[cur_pos:]))
	elif list_analize_type=="Evkl":
		result.append([np.linalg.norm(np.array([a/s1 for a in array1[cur_pos:]])-np.array([b/s2 for b in array2[cur_pos:]])),0])
	elif list_analize_type=="proc":
		result.append(compare_nonzerro(array1[cur_pos:],array2[cur_pos:]))
	else:
		raise "Unknown list_analize_type"


  return result
  
def do_correction_f(data,reference): #This function corrects data according to standart error of reference
	reference_y=[t[0] for t in reference]
	data_y=[t[0] for t in data]

#	ma_ref=max([e for e in e1 for e1 in reference_y])
	ma_ref=1.0 # Let's say correlation for reference schould be 1 everywhere
	mi_ref=min([float(e) for e in reference_y if e != 'Undef'])

	ma_dat=max([float(e) for e in data_y if e != 'Undef'])
	mi_dat=min([float(e) for e in data_y if e != 'Undef'])

	sd_ref=np.std(np.array([float(e) for e in reference_y if e != 'Undef']))
	mean_ref=np.mean(np.array([float(e) for e in reference_y if e != 'Undef']))
	result=[]
	coefficient = (float(ma_dat)-mi_dat)/(ma_ref-mi_ref)	

	print "borders for correction = ", mean_ref, " +/- ", 1.0*sd_ref, "min= ",mi_ref,"max= ",ma_ref
	print "new_min= ",min([float(e) for ind,e in enumerate(reference_y) if (e != 'Undef') and (data_y[ind] != 'Undef')])
	removed_dots=0
	ind_i=0
	while ind_i<len(data_y):
#		if (data_y[ind_i] == 'Undef') or (reference_y[ind_i] == 'Undef'):
#			if len(result) > 3:
#				result[-1] = ['Undef','Undef']
#				result[-2] = ['Undef','Undef']
#				result[-3] = ['Undef','Undef']
#			result.append(['Undef','Undef'])
#			result.append(['Undef','Undef'])
#			result.append(['Undef','Undef'])
#			ind_i += 2
#		else:
		if (reference_y[ind_i] == 'Undef') or abs(reference_y[ind_i]-mean_ref) > 1.0*sd_ref:
				for temp_index in xrange(max(len(result)-1,0),len(result)):
					result[temp_index] = ['Undef','Undef']
					removed_dots += 1
				for temp_index in xrange(ind_i,min(ind_i+1,len(data_y)-1)):
					result.append(['Undef','Undef'])
					removed_dots += 1
					ind_i += 1
		else:
				result.append([data_y[ind_i],data[ind_i]])
				ind_i += 1
	#			result.append([float(data_y[ind_i])+(ma_ref-float(reference_y[ind_i]))*coefficient,max(data[ind_i],reference[ind_i])])
	print removed_dots, " dots removed during correction"
	return result

def write_cor_to_file(di_1,di_2,o_filename,chrms_in_bins,chrms_axes,genome_db,
			window_size=2,clear_plot=True,
			annot="Correlation",
			color=0,
			list_analize_type="Sp",
			do_correction=False,correctionlist=[]): #do correction - coorect all values according to some graph

  if hasattr(di_1[0], '__iter__'):
	comment=""
  else:
	comment="Whole list correlation = "+str(scipy.stats.spearmanr(di_1,di_2))
  
  comment += " Window="+str(window_size)
  
  filename_extenson = ""
  if do_correction: filename_extenson += ".corrected"
  filename_extenson += ".bed"
  
  f3=open(o_filename+str(window_size)+annot+"."+list_analize_type+filename_extenson,'w')
  header="track name='"+annot+"."+list_analize_type+"' description='"+comment+"' visibility=full color=200,100,0 altColor=0,100,200 priority=20\n"
  f3.write(header)

  t=window_sperman(di_1,di_2,int(window_size),list_analize_type)
  if do_correction:
	t=do_correction_f(t,correctionlist)

  
  count=0
  
  plot_SpCorValues_X_array=[]
  plot_SpCorValues_X_array.append([])
  plot_Chromosomes_X_array=[]
  plot_chr_Y_array={} # Y-values for chromosomes line
  plot_SpCorValues_Y_array=[]
  plot_SpCorValues_Y_array.append([])
  for i in t:
	if i[0] != 'Undef': # Just skip zerro-bins in graph #!!!!!!!!!!!!!!!!!!!-Change it later-!!!!!!!!!!!!!
		plot_SpCorValues_Y_array[-1].append(float(i[0]))
		plot_SpCorValues_X_array[-1].append(count+1)
	elif len(plot_SpCorValues_Y_array[-1])>0: #split line into fragments, where values are defined (no zero bins)
		plot_SpCorValues_Y_array.append([])
		plot_SpCorValues_X_array.append([])
   
	#This part makes chromosomes ready to be drawn later
	curChrmIdx=chrms_in_bins[int(window_size)*int(count)]
	if (not curChrmIdx in plot_chr_Y_array.keys()):
		plot_chr_Y_array[curChrmIdx]=[]
	plot_chr_Y_array[curChrmIdx].append(count+1) # Y-values for chromosomes line
	plot_Chromosomes_X_array.append(count+1)
   
	#Save data to file
	if curChrmIdx < 19: 
		chr_to_write='chr'+str(curChrmIdx+1) 
	else: 
		chr_to_write='chrX'
	if curChrmIdx == 0: 
		curRelativeBinNumb = count*int(window_size) 
	else: 
		curRelativeBinNumb = count*int(window_size)-genome_db.chrmLensBin[0:curChrmIdx].sum()

	strToWrite = str(chr_to_write)+"\t" #cur chromosome in bed numeration format
	strToWrite += str(genome_db.posBinCont[count*int(window_size)])+"\t" #start nucleotide of current bin
	strToWrite += str(genome_db.posBinCont[count*int(window_size)]+genome_db.binSizesBp[curChrmIdx][curRelativeBinNumb]) #end nucleotide of current bin
	if i[0] != 'Undef': f3.write(strToWrite+"\t"+str(i[0])+"\n")
	
	#iterate count
	count=count+1
   
  #check that last fragment is not empty
  if len(plot_SpCorValues_Y_array[-1])==0:
	del plot_SpCorValues_Y_array[-1]
	del plot_SpCorValues_X_array[-1]

  # max Y-value for chromosomes graph is minimum value for data graph
  Chr_plot_maxVal=min([min(m) for m in plot_SpCorValues_Y_array])
    
  
  if (clear_plot):
   plt.clf()  #clear image
   
  if (not chrms_axes):

   #format axes
   Chr_plot_minVal=-1
   
   ylabel=""
   if len(list_analize_type.split("_w"))>0:
	list_analize_type=list_analize_type.split("_w")[0]
#	ylabel=" (inner_window_size="+str(inner_window_size)+")"
   if list_analize_type=="Sp":
	plt.ylabel('Sperman correlation coefficient'+ylabel)
   elif list_analize_type=="pears":
	plt.ylabel('Pearson correlation coefficient'+ylabel)
   elif list_analize_type=="Evkl":
	plt.ylabel('Euclidean distance'+ylabel)
   elif list_analize_type=="proc":
    plt.ylabel('Manhattan distance'+ylabel)
   
   if window_size == 1:
	plt.xlabel('Window number')
   else:
    plt.xlabel('Bin number')   
   max_y_val=max([max(m) for m in plot_SpCorValues_Y_array])
   max_y_val += float(max_y_val)/10
   plt.axis([0,count+2,Chr_plot_minVal-0.2,max_y_val ])
   ax=plt.gca()
   ax.set_yticks(np.arange(Chr_plot_minVal-0.2,max_y_val,0.1))
   if count < 101 :
	ax.set_xticks(np.arange(0,count))
   else:
	ax.set_xticks(np.arange(0,count,round(count/100)*10))
  
   ax.yaxis.grid(b=True, which='both', color='gray', linestyle='dotted',linewidth=0.5)
   ax.xaxis.grid(b=False)
   
   #Print chromosomes in the bottom of graph
   
   for i in plot_chr_Y_array.keys(): 
	x_range=[]
	y_range=[]

	y_val=Chr_plot_minVal+float(i)*abs(Chr_plot_maxVal-Chr_plot_minVal)/len(plot_chr_Y_array.keys())
	xmin_val=min(plot_chr_Y_array[i])
	xmax_val=max(plot_chr_Y_array[i])
	x_range=range(xmin_val,xmax_val+1)
	for j in x_range:
		y_range.append(y_val)
	plt.plot(x_range,y_range)
	plt.vlines(xmax_val+1, Chr_plot_minVal-0.5, 1, linestyles='dotted',linewidth=0.5)
		
  line_format=colors[color]+markers[color]+'-' #format for main line
  for i in range(len(plot_SpCorValues_Y_array)-1):
#	print '>>>>>>>>>>>>>>',i,"\n",plot_SpCorValues_X_array[i],plot_SpCorValues_Y_array[i]
	plt.plot(plot_SpCorValues_X_array[i],plot_SpCorValues_Y_array[i],line_format,markersize=1,color=colors[color])
  plt.plot(plot_SpCorValues_X_array[-1],plot_SpCorValues_Y_array[-1],line_format,markersize=1,color=colors[color],label=annot)
  plt.legend(loc=4)

  f3.close()
  return t
