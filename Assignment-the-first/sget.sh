#!/bin/bash 


#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --output=index2_hist.%j.out
#SBATCH --cpus-per-task=8

conda activate bgmp_py310



# /usr/bin/time -v get_hist.py \


/usr/bin/time -v ./get_hist.py \
    -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
    -o index2 \
    -l 8




