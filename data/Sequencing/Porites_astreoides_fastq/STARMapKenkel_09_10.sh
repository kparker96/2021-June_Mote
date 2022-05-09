#!/bin/bash -l

#SBATCH -o 2022-04-21_STARMapKenkel_09_10.txt
#SBATCH -n 4
#SBATCH --mail-user=kpark049@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=STARKenkel_09_10

enable_lmod

module load container_env star

STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_30_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_30_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_30_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_30_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_30_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_30_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_30_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_30_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_30_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_30_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_30_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_30_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_30_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_30_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_30_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_34_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_34_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_34_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_34_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_34_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_34_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_34_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_34_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_34_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_34_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_34_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_34_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_34_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_34_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_34_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_37_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_37_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_37_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_37_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_37_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_37_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_37_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_37_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_39_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_39_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_39_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_39_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_39_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_39_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_39_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_39_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_39_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_39_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_39_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_39_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_39_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_39_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_39_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_only_seq_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_only_seq_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_only_seq_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_T0_0_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_T0_0_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_T0_0_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_34_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_34_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_34_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_34_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_37_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_37_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_37_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_39_1140_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_39_180_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_39_360_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_39_750_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_only_seq_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_only_seq_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_only_seq_RNASeq_2kenkel_
STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_2kenkel_
