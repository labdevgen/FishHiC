import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)

if not os.path.exists('../../data/serov/tmp/'):
    os.mkdir('../../data/serov/tmp/')

# A. Map the reads iteratively.
mapping.iterative_mapping(
    bowtie_path='../../bin/bowtie2/bowtie2',
    bowtie_index_path='../../bin/bowtie2/index/mm10',
    fastq_path='../../data/serov/HiC_Fib.fastq',
    out_sam_path='../../data/serov/HiC_Fib_1.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=0,
    seq_end=50,
    nthreads=4,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
    temp_dir='../../data/serov/tmp',  # optional, keep temporary files here
    bowtie_flags='--very-sensitive')
 #   bash_reader='../../bin/sra/bin/fastq-dump -Z')

mapping.iterative_mapping(
    bowtie_path='../../bin/bowtie2/bowtie2',
    bowtie_index_path='../../bin/bowtie2/index/mm10',
    fastq_path='../../data/serov/HiC_Fib.fastq',
    out_sam_path='../../data/serov/HiC_Fib_2.bam',
    min_seq_len=25,
    len_step=5,
    seq_start=51,
    seq_end=100,
    nthreads=4,  
    #max_reads_per_chunk = 10000000, 
    temp_dir='../../data/serov/tmp',  
    bowtie_flags='--very-sensitive')
   # bash_reader='../../bin/sra/bin/fastq-dump -Z')

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
mapped_reads_Fib = h5dict.h5dict('../../data/serov/mapped_reads_Fib.hdf5')
genome_db    = genome.Genome('../../fasta/mm10', readChrms=['#', 'X'])

mapping.parse_sam(
    sam_basename1='../../data/serov/HiC_Fib_1.bam',
    sam_basename2='../../data/serov/HiC_Fib_2.bam',
    out_dict=mapped_reads_Fib,
    genome_db=genome_db, 
    enzyme_name='HindIII')

