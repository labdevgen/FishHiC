import os
import logging
from hiclib import mapping
from mirnylib import h5dict, genome

logging.basicConfig(level=logging.DEBUG)

# B. Parse the mapped sequences into a Python data structure,
#    assign the ultra-sonic fragments to restriction fragments.
mapped_reads_Sp1 = h5dict.h5dict('../../data/serov/mapped_reads_Sp1.hdf5')
genome_db    = genome.Genome('../../fasta/mm10', readChrms=['#', 'X'])

mapping.parse_sam(
    sam_basename1='../../data/serov/HiC_Sp1_1.bam',
    sam_basename2='../../data/serov/HiC_Sp1_2.bam',
    out_dict=mapped_reads_Sp1,
    genome_db=genome_db, 
    enzyme_name='HindIII')

