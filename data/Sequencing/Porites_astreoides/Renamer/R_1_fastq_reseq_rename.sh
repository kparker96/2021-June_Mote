#!/bin/bash -l

#SBATCH -o reseq_fastq_rename_R1.txt
#SBATCH -n 1
#SBATCH --mail-user=kpark049@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=reseqFastqR1

/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_1_fastq_reseq_rename.txt