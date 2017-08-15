import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)

#FASTQ_fpath='/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/SRR443883.sra'
FASTQ_fpath='/ifs/DATA/opistorchis/Fishman/data/Serov/130422_GA470.HiC_Fib.read0.fastq'

#out_sam_fpath='/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/SRR443883'
out_sam_fpath='/ifs/DATA/opistorchis/Fishman/data/Serov/130422_GA470.HiC_Fib.read0'

temp_dir_path='/ifs/DATA/opistorchis/Fishman/tmp/'
genome_name='mm10'

fastq_path='/ifs/home/bionet/pomaznoy/soft/HiC/mirnylab-hiclib-04615fd1aece/bin/sra/bin/'

if not os.path.exists(temp_dir_path):
    os.mkdir(temp_dir_path)

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='../../bin/bowtie2/bowtie2',
    bowtie_index_path='../../bin/bowtie2/index/'+genome_name,
    fastq_path=FASTQ_fpath,
    out_sam_path=out_sam_fpath+'_1.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=0,
    seq_end=49,
    nthreads=8,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir=temp_dir_path,  # optional, keep temporary files here
    bowtie_flags='--very-sensitive',
    bash_reader='cat')  #fastq_path+'fastq-dump -Z')

mapping.iterative_mapping(
    bowtie_path='../../bin/bowtie2/bowtie2',
    bowtie_index_path='../../bin/bowtie2/index/'+genome_name,
    fastq_path=FASTQ_fpath,
    out_sam_path=out_sam_fpath+'_2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=50,
    seq_end=99,
    nthreads=8,  
    #max_reads_per_chunk = 10000000, 
    temp_dir=temp_dir_path,  
    bowtie_flags='--very-sensitive',
    bash_reader='cat') #fastq_path+'fastq-dump -Z')

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
mapped_reads = h5dict.h5dict('/ifs/DATA/opistorchis/Fishman/data/Sample/mESC/mapped_reads.hdf5')
genome_db    = genome.Genome('../../fasta/'+genome_name, readChrms=['#', 'X'])

mapping.parse_sam(
    sam_basename1=out_sam_fpath+'_1.bam',
    sam_basename2=out_sam_fpath+'_2.bam',
    out_dict=mapped_reads,
    genome_db=genome_db, 
    enzyme_name='HindIII')

