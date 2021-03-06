import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)

base_folder='/mnt/storage/home/vsfishman/HiC/data/'
base_filename = "Fib_full2"
tmp_folder='/mnt/storage/home/vsfishman/tmp/HiC_tmp'
maped_reads_filepath=base_folder+'mapped_reads_'+base_filename+'.hdf5'

FASTQ_fpath="/mnt/storage/home/vsfishman/tmp/Distr/Fib_full2.fastq"
out_sam_fpath=tmp_folder+'/'+base_filename
genome_name='mm9'

if not os.path.exists(tmp_folder):
    os.mkdir(tmp_folder)

# A. Map the reads iteratively.
#mapping.iterative_mapping(
#    bowtie_path='../bin/bowtie2/bowtie2',
#    bowtie_index_path='../bin/bowtie2/index/'+genome_name,
#    fastq_path=FASTQ_fpath,
#    out_sam_path=out_sam_fpath+'_1.bam',
#    min_seq_len=25,
#    len_step=5,
#    seq_start=0,
#    seq_end=49,
#    nthreads=8,  # on intel corei7 CPUs 4 threads are as fast as
                 # 8, but leave some room for you other applications
    #max_reads_per_chunk = 10000000,  #optional, on low-memory machines
#    temp_dir=tmp_folder,  # optional, keep temporary files here
#    bowtie_flags='--very-sensitive',
#   bash_reader=None)#../../bin/sra/bin/fastq-dump -Z')

#mapping.iterative_mapping(
#    bowtie_path='../bin/bowtie2/bowtie2',
#    bowtie_index_path='../bin/bowtie2/index/'+genome_name,
#    fastq_path=FASTQ_fpath,
#    out_sam_path=out_sam_fpath+'_2.bam',
#    min_seq_len=25,
#    len_step=5,
#    seq_start=50,
#    seq_end=99,
#    nthreads=8,  
#    #max_reads_per_chunk = 10000000, 
#    temp_dir=tmp_folder,  
#    bowtie_flags='--very-sensitive',
#    bash_reader=None)#../../bin/sra/bin/fastq-dump -Z')

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
mapped_reads = h5dict.h5dict(maped_reads_filepath)
genome_db    = genome.Genome('../fasta/'+genome_name, readChrms=['#', 'X'])

mapping.parse_sam(
    sam_basename1=out_sam_fpath+'_1.bam',
    sam_basename2=out_sam_fpath+'_2.bam',
    out_dict=mapped_reads,
    genome_db=genome_db, 
    enzyme_name='HindIII', save_seqs=True)
