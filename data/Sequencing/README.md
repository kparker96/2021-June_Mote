# *Porites astreoides* RNAseq Data Analysis Log
---
## 2022-01-30: Created Repository and Imported File Name List Locally
 Raw fq.gz files were moved to this directory:  
	
	$pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

A list of raw files names was made for forward and reverse reads:  
	- raw\_fastqs\_names\_1.txt  
	- raw\_fastqs\_names\_2.txt

	$head raw_fastqs_names_1.txt 
	R100_1.fq.gz
	R101_1.fq.gz
	R10_1.fq.gz
	R102_1.fq.gz
	R103_1.fq.gz
	R104_1.fq.gz
	R105_1.fq.gz
	R106_1.fq.gz
	R107_1.fq.gz
	R108_1.fq.gz

	$head raw_fastqs_names_2.txt
	R100_2.fq.gz
	R101_2.fq.gz
	R102_2.fq.gz
	R10_2.fq.gz
	R103_2.fq.gz
	R104_2.fq.gz
	R105_2.fq.gz
	R106_2.fq.gz
	R107_2.fq.gz
	R108_2.fq.gz

Both files were copied locally to:  
C:\Users\kpark\OneDrive\Documents\Barshis_Lab\2021-June_Mote\data\Sequencing\Raw_Fastqs
	
	$scp kpark049@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/raw_fastqs_names_*.txt ./

Using the Notepad ++ "find and replace" function, "\_1.fq.gz" or "\_1.fq.gz" was removed from each file and the header "fastq" was added. These new files saved locally and named:  
	- raw\_fastqs\_sample\_names\_1.txt  
	- raw\_fastqs\_sample\_names\_2.txt

	Example:
	fastq
	R100
	R101
	R10
	R102
	...

## 2022-01-31: Sample Name Table Made to Run renamer.py

The fastq sample name lists 1 and 2 were imported into R (2021-June-Mote.Rproj) along with a table of the sample names with corresponding sample ID. 

	# Raw list from cluster after find and replacing "_1.fq.gz"
	R_1 <- read.delim("data/Sequencing/Raw_Fastqs/raw_fastqs_sample_names_1.txt")

	# Raw list from cluster after find and replacing "_2.fq.gz"
	R_2 <- read.delim("data/Sequencing/Raw_Fastqs/raw_fastqs_sample_names_2.txt")

	# Sample_ID list taken from datasheet returned from novogene 
	sample_names <- read.delim("data/Sequencing/Raw_Fastqs/sample_names.txt")

	setdiff(R_1$fastq, sample_names$fastq) 

sample\_names.txt
  
	> head(sample_names)
	Sample_ID fastq
	1 IPa06-C-TF    R1
	2 IPa08-C-TF    R3
	3 IPa09-C-TF    R4
	4 IPa07-C-T0    R7
	5 IPa08-C-T0    R8
	6 IPa06-C-T1   R11

The fastq sample lists were combined with the sample names table to make sure the fastq sample names lined up with the correct sample ID. 

	# R_1 TABLE #  

	# Joining the raw cluster list with the sample IDs for R_1 and adding _1.fq.gz file name
	R_1_rename <- right_join(R_1, sample_names, by = "fastq") %>%
	mutate(old_fastq = paste(fastq, "_1.fq.gz", sep = "")) %>%
	mutate(new_fastq = paste(Sample_ID, "_1.fq.gz", sep = ""))

	# Final Sample Table for R_1
	R_1_rename_ready <- R_1_rename %>%
	select(old_fastq, new_fastq)

	# Export Table 
	write.table(R_1_rename_ready, file = "data/Sequencing/Raw_Fastqs/R_1_fastq_rename.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


	# R_2 TABLE # 

	# Joining the raw cluster list with the sample IDs for R_1 and adding _2.fq.gz file name
	R_2_rename <- right_join(R_2, sample_names, by = "fastq") %>%
	mutate(old_fastq = paste(fastq, "_2.fq.gz", sep = "")) %>%
	mutate(new_fastq = paste(Sample_ID, "_2.fq.gz", sep = ""))

	# Final Sample Table for R_2
	R_2_rename_ready <- R_2_rename %>%
	select(old_fastq, new_fastq)

	# Export Table 
	write.table(R_2_rename_ready, file = "data/Sequencing/Raw_Fastqs/R_2_fastq_rename.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

	
Exported tables were then copied into the directory on the cluster.  

**NOTE:** Make sure to add headers in the future. I didn't do that here and the first line file name did not work with the renamer.py because it was considered the header. 

	$pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$head R_1_fastq_rename.txt 
	R100_1.fq.gz    IPa10-L-T5_1.fq.gz
	R101_1.fq.gz    IPa06-M-T5_1.fq.gz
	R10_1.fq.gz     IPa10-C-T0_1.fq.gz
	R102_1.fq.gz    IPa07-M-T5_1.fq.gz
	R103_1.fq.gz    IPa08-M-T5_1.fq.gz
	R104_1.fq.gz    IPa09-M-T5_1.fq.gz
	R105_1.fq.gz    IPa10-M-T5_1.fq.gz
	R106_1.fq.gz    IPa06-H-T5_1.fq.gz
	R107_1.fq.gz    IPa07-H-T5_1.fq.gz
	R108_1.fq.gz    IPa08-H-T5_1.fq.gz

	$head R_2_fastq_rename.txt
	R100_2.fq.gz    IPa10-L-T5_2.fq.gz
	R101_2.fq.gz    IPa06-M-T5_2.fq.gz
	R102_2.fq.gz    IPa07-M-T5_2.fq.gz
	R10_2.fq.gz     IPa10-C-T0_2.fq.gz
	R103_2.fq.gz    IPa08-M-T5_2.fq.gz
	R104_2.fq.gz    IPa09-M-T5_2.fq.gz
	R105_2.fq.gz    IPa10-M-T5_2.fq.gz
	R106_2.fq.gz    IPa06-H-T5_2.fq.gz
	R107_2.fq.gz    IPa07-H-T5_2.fq.gz
	R108_2.fq.gz    IPa08-H-T5_2.fq.gz


A renamer script was made in my scripts folder.  
**NOTE:** un-comment the "print" line and comment the "cols" line to get a preview of the rename. 

	$pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts

	$ nano renamer_KEP.py
	#!/usr/bin/env python
	####usage renamer.py renamingtable
	#### this script take the entries in the first column of table and renames (mv's) them to files with the $import sys
	import os

	fin=open(sys.argv[1],'r')
	linecount=0
	for line in fin:
        	linecount+=1
        	if linecount>=2:
               	cols=line.rstrip().split('\t')
      #			print 'mv %s %s' %(cols[0], cols[1])
                os.popen('mv %s %s' %(cols[0], cols[1]))

Sbatch scripts were made in the working directory to run the renmaer_KEP.py script.

	$pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/
	 
	$nano KparkerRename_R_1.sh
	#!/bin/bash -l


	#SBATCH -o raw_fastq_rename_R1.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=KPrenameFastq


	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_1_rextract_names.txt


	$nano KparkerRename_R_2.sh
	#!/bin/bash -l


	#SBATCH -o raw_fastq_rename_R1.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=KPrenameFastq


	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_2_rextract_names.txt
	
	
Sbatch scripts were then run for forward and reverse raw .fq.gz files. 

	$ salloc
	salloc: Pending job allocation 9693629
	salloc: job 9693629 queued and waiting for resources
	salloc: job 9693629 has been allocated resources
	salloc: Granted job allocation 9693629

	$sbatch KparkerRename_R_1.sh
	

	$sbatch KparkerRename_R_2.sh


After the job finished most of the files were renamed. 

	Example: 
	
	$ls

	2MAc48-M-T5_1.fq.gz   IPa08-C-T5_2.fq.gz  OPa-01-C-TF_2.fq      OPa-04-H-T4_2.fq.gz
	2MAc48-M-T5_2.fq.gz   IPa08-C-TF_1.fq.gz  OPa-01-H-T1_1.fq      OPa-04-H-T5_1.fq.gz
	2MAc50-H-T2_1.fq.gz   IPa08-C-TF_2.fq.gz  OPa-01-H-T1_2.fq      OPa-04-H-T5_2.fq.gz
	2MAc50-H-T2_2.fq.gz   IPa08-H-T1_1.fq.gz  OPa-01-H-T2_1.fq      OPa-04-L-T1_1.fq.gz
	2MAc50-H-T3_1.fq.gz   IPa08-H-T1_2.fq.gz  OPa-01-H-T2_2.fq      OPa-04-L-T1_2.fq.gz
	2MAc50-H-T3_2.fq.gz   IPa08-H-T2_1.fq.gz  OPa-01-H-T3_1.fq      OPa-04-L-T2_1.fq.gz
	...

~10 samples did not get renamed because the tables used did not have headers and the lowercase r was missing in the table for re-extracted samples that were sequenced. New sample name tables were made for the missed samples with headers and the sbatch script was run again. 

## 2022-02-01