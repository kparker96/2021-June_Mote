#!/bin/bash -l

#SBATCH -o reseq_fastq_rename_R2.txt
#SBATCH -n 1
#SBATCH --mail-user=kpark049@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=reseqFastqR2

/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_2_fastq_reseq_rename.txt