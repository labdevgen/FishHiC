import datetime
import random

now = datetime.datetime.now()

print "Starting domains_search_Dixon "
print str(now)
print "\n"

import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

from mirnylib import genome
from mirnylib import h5dict
from mirnylib import plotting
from hiclib import binnedData
from hiclib import fragmentHiC

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--hmap")
args = parser.parse_args()

heatmap_filepath=str(args.hmap).strip()

def get_domain_res(hmap):
  raw_heatmap = h5dict.h5dict(hmap, mode='r') 
  domain_res = int(raw_heatmap['resolution'])
  print "resolution defined by heatmap: ",domain_res
  return domain_res

domain_res = get_domain_res(heatmap_filepath)

print "Domain resolution = "+str(domain_res)
base_filename = heatmap_filepath.split("/")[-1].replace("-","_")


#Define files and directories
base_folder='/mnt/storage/home/vsfishman/HiC/data/chick/DixonDomains'
HMM_file_path='/mnt/storage/home/vsfishman/HiC/DI/domaincall_software/'
DI_from_matrix_file_path = "/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/DI_from_matrix_minja.pl"


genome_db = genome.Genome("/mnt/storage/home/vsfishman/HiC/fasta/GalGal5/GCF_000002315.4_Gallus_gallus-5.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/",
				readChrms=[],
				chrmFileTemplate="%s.fna")

genome_fai_filepath=genome_db.genomePath+'/GalGal5ChrmLevel.fai'

def executeBashCommand(bashCommand,doPrint=True):
  print "executing command %s \n" % (bashCommand)
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  if doPrint:
    print p.communicate()[0]
    return ""
  else:
    return p.communicate()[0]

def doDixonSearch():
  print "----Starting dixon domain search---\n"
  bashCommand =DI_from_matrix_file_path+" "+heatmap_filepath+".matrix "+str(domain_res)+" "+str(domain_res*50)+" "+genome_fai_filepath+" "+heatmap_filepath+".matrix.di"
  print "Executing ",bashCommand
  print executeBashCommand(bashCommand)
  print "\n Done!\n"

def splitChrmsColFile():
  ChrmsCol_file_name=heatmap_filepath+".matrix.di.m.hmm_7colfile"
  executeBashCommand ("mkdir -p "+ base_folder + base_filename + "_domains_"+str(domain_res/1000)+'KB') 
  print "Opening ChrmsCol_file " + ChrmsCol_file_name+"\n"

  f_in=open(ChrmsCol_file_name,'r')
  chr_numb="chr0"
  f_out=open(base_folder + base_filename + "_domains_"+str(domain_res/1000)+'KB/'+base_filename+"_"+str(chr_numb)+".col7_out",'w')
  chrms_list=[chr_numb]
  for line in f_in:
    new_chr_numb=line.split("\t")[0]
    if new_chr_numb != chr_numb:
       chr_numb=new_chr_numb
       f_out.close()
       chrms_list.append(chr_numb)
       f_out=open(base_folder + base_filename + "_domains_"+str(domain_res/1000)+'KB/'+base_filename+"_"+str(chr_numb)+".col7_out",'w')
    f_out.write(line)

  for chr_numb in chrms_list:
#        print "old command="+"perl "+HMM_file_path+"perl_scripts/hmm_probablity_correcter.pl "+base_folder+base_filename+"_domains/"+base_filename+str(chr_numb)+".col7_out "+str(2)+" "+str(0.99)+" "+str(domain_res)+" | perl "+HMM_file_path+"perl_scripts/hmm-state_caller.pl "+genome_fai_filepath+" "+chr_numb.split("chr")[1]+" | perl "+HMM_file_path+"perl_scripts/hmm-state_domains.pl > "+base_folder+base_filename+"_domains/"+base_filename+chr_numb+".final_domains\n"

	bashCommand="perl "+HMM_file_path+"perl_scripts/hmm_probablity_correcter.pl "+base_folder + base_filename + "_domains_"+str(domain_res/1000)+'KB/'+base_filename+"_"+str(chr_numb)+".col7_out "+str(2)+" "+str(0.99)+" "+str(domain_res)

#+" | perl "+HMM_file_path+"perl_scripts/hmm-state_caller.pl "+genome_fai_filepath+" "+chr_numb.split("chr")[1]
#+" | perl "+HMM_file_path+"perl_scripts/hmm-state_domains.pl > "+base_folder+base_filename+"_domains/"+base_filename+chr_numb+".final_domains"

	p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	t_out=p.communicate()[0]

	bashCommand="perl "+HMM_file_path+"perl_scripts/hmm-state_caller.pl "+genome_fai_filepath+" "+chr_numb.split("chr")[1]
	p=subprocess.Popen(bashCommand.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	t_out=p.communicate(input=t_out)[0]
	
	bashCommand="perl "+HMM_file_path+"perl_scripts/hmm-state_domains.pl"
	p=subprocess.Popen(bashCommand.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	t_out=p.communicate(input=t_out)[0]
	f_out=open(base_folder + base_filename + "_domains_"+str(domain_res/1000)+'KB/'+base_filename+"_"+str(chr_numb)+".final_domains",'w')
	f_out.write(t_out)
	f_out.close()	
  print "\n Done!\n";

def doDixon_PostProc():
  print "----Starting dixon postproc---\n"
  if (not os.path.isfile(heatmap_filepath+".matrix.di.m")):
    print "Generating new HMM_calls file ("+HMM_file_path+'HMM_calls_'+base_filename+".m)\n"
    m_inp_file = open('HMM_calls.m','r')
    fname = 'HMM_calls_'+str(random.randint(0,1000))
    m_out_file = open(fname+".m",'w')
    m_content = m_inp_file.read()
    m_content = m_content.replace('i_filename',heatmap_filepath+".matrix.di")
    m_content = m_content.replace("addpath(genpath('","addpath(genpath('"+HMM_file_path)

    m_out_file.write(m_content)
    m_inp_file.close()
    m_out_file.close()
#    executeBashCommand("bash ~/modules.sh'")
    executeBashCommand("bash -c 'module load matlab/r2013b'")
    bashCommand ="matlab nodisplay -nodesktop -nojvm -nosplash -r "+fname
    executeBashCommand(bashCommand)
#    executeBashCommand("rm "+fname)

  print ".m file is ready \n"
  bashCommand="perl "+HMM_file_path+"perl_scripts/file_ends_cleaner.pl "+heatmap_filepath+".matrix.di.m "+heatmap_filepath+".matrix.di"  
  p=subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  t_out=p.communicate()[0]

  bashCommand="perl "+HMM_file_path+"perl_scripts/converter_7col.pl"
  p=subprocess.Popen(bashCommand.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
  t_out=p.communicate(input=t_out)[0]

  f_out=open(heatmap_filepath+".matrix.di.m.hmm_7colfile",'w')
  f_out.write(t_out)
  f_out.close()	
  print "\n Done!\n"
  

def CreateMatrixFile():

  BD = binnedData.binnedData(domain_res, genome_db)
  BD.simpleLoad(heatmap_filepath, 'heatmap')
  print "Writing file %s with the information \n"% (heatmap_filepath+'.matrix')
  print "Format:\nChrIndex \t StartBin(Nucleotyde) \t EndBin(Nucl) \t Values\n"
  f=open(heatmap_filepath+'.matrix','w')
  for i in range (len(BD.dataDict['heatmap'])):
    strToWrite = ""
    curChrmIdx = genome_db.chrmIdxBinCont[i]
    if curChrmIdx == 0: 
      curRelativeBinNumb = i 
    else: 
      curRelativeBinNumb = i-genome_db.chrmLensBin[0:curChrmIdx].sum()

    strToWrite += str(curChrmIdx)+"\t"+str(genome_db.posBinCont[i])+"\t"+str(genome_db.posBinCont[i]+genome_db.binSizesBp[curChrmIdx][curRelativeBinNumb])
    for j in range (len(BD.dataDict['heatmap'])): 
       strToWrite += "\t"+str(BD.dataDict['heatmap'][i][j])
    strToWrite += "\n"
    f.write (strToWrite)
  f.close()


if (not os.path.isfile(heatmap_filepath+".matrix")):
   print "Generating matrix file %s\n" % (heatmap_filepath+".matrix")
   CreateMatrixFile()
   doDixonSearch()
   doDixon_PostProc()
   splitChrmsColFile()
else:
   print "Found file ",heatmap_filepath+".matrix"
   
if (not os.path.isfile(heatmap_filepath+".matrix.di")):
   doDixonSearch()
   doDixon_PostProc()
   splitChrmsColFile()
else:
   print "Found file ",heatmap_filepath+".matrix.di"
 
if (not os.path.isfile(heatmap_filepath+".matrix.di.m.hmm_7colfile")):
   doDixon_PostProc()
   splitChrmsColFile()

splitChrmsColFile()
print "Done!\n"
