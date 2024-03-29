#!/bin/bash -l

#SBATCH -o 2022-04-14_STARMapKenkel_01.txt
#SBATCH -n 4
#SBATCH --mail-user=kpark049@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=STARKenkel_01

enable_lmod

module load container_env star

for i in 20210608T1400_CBASS_US_MoteOff_Past_1_*_RNASeq_R1_val_1.fq.gz ; do `STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal  --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn $i ${i%_R1_val_1.fq.gz}_R2_val_2.fq.gz --outFileNamePrefix ${i%_R1_val_1.fq.gz}_2kenkel`; done
