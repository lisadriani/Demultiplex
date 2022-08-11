#!/bin/bash 


#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --output=demultiplex.%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=ladriani@uoregon.edu  ### Where to send mail
#SBATCH --job-name=demultiplex         ### Job Name


conda activate bgmp_py310

/usr/bin/time -v ./the_code.py \
    -R1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
    -R2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz\
    -I1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
    -I2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
    -I indexes.txt


