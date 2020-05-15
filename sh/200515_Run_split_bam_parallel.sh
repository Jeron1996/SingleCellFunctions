#!/bin/bash -e

#This submits the split bam script to a SGE scheduler
qsub -cwd -V -N split_parallel -pe smp 10 -l mem_requested=10G,tmp_requested=50G -b y -j y './200515_split_bam_parallel.sh'
