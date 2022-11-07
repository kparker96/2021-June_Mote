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

	$nano renamer_KEP.py
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

## 2022-02-01: Restarting the Pipeline, Renamer.py

All previous work was moved to the new directory "old". The raw zip file of raw data was re-copied into my directory. 
	
	$pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1
	
	$ mkdir old

	$mv raw_data_fastqs/ X204SC21081158-Z01-F002/ old

	$ls
	old  X204SC21081158-Z01-F002.zip

	$salloc
	salloc: Pending job allocation 9693669
	salloc: job 9693669 queued and waiting for resources
	salloc: job 9693669 has been allocated resources
	salloc: Granted job allocation 9693669
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.

	$pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1

	$nano unzip_raw.sh
	#!/bin/bash -l

	#SBATCH -o KPunzip.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=KPunzip

	unzip X204SC21081158-Z01-F002.zip

	$sbatch unzip_raw.sh
	Submitted batch job 9693670

A new "raw\_data\_fastqs"directory was made. 

	$pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1
	
	$mkdir raw_data_fastqs
	
	$ls 
	KPunzip.txt  __MACOSX  old  raw_data_fastqs  unzip_raw.sh  X204SC21081158-Z01-F002  X204SC21081158-Z01-F002.zip

All raw data files from X204SC21081158-Z01-F002/ were moved into the raw\_data\_fastqs directory. 

	$pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/

	$mv /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/*.fq.gz ./

All files were renamed to informative sample ID naming. 

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$ salloc
	salloc: Pending job allocation 9693775
	salloc: job 9693775 queued and waiting for resources
	salloc: job 9693775 has been allocated resources
	salloc: Granted job allocation 9693775
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.

	$ cat KparkerRenamer_R_1.sh
	#!/bin/bash -l

	#SBATCH -o raw_fastq_rename_R1.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=KPrenameFastq

	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_1_fastq_rename.txt

	$ head R_1_fastq_rename.txt
	OldName NewName
	R100_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1.fq.gz
	R101_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R1.fq.gz
	R10_1.fq.gz     20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R1.fq.gz
	R102_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R1.fq.gz
	R103_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R1.fq.gz
	R104_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R1.fq.gz
	R105_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R1.fq.gz
	R106_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_6_39_1140_RNASeq_R1.fq.gz
	R107_1.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_7_39_1140_RNASeq_R1.fq.gz

	$ sbatch KparkerRenamer_R_1.sh
	Submitted batch job 9693776

	$ cat KparkerRenamer_R_2.sh
	#!/bin/bash -l

	#SBATCH -o raw_fastq_rename_R1.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=KPrenameFastq

	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_2_fastq_rename.txt

	$ head R_2_fastq_rename.txt
	OldName NewName
	R100_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2.fq.gz
	R101_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R2.fq.gz
	R10_2.fq.gz     20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R2.fq.gz
	R102_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R2.fq.gz
	R103_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R2.fq.gz
	R104_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R2.fq.gz
	R105_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R2.fq.gz
	R106_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_6_39_1140_RNASeq_R2.fq.gz
	R107_2.fq.gz    20210608T1400_CBASS_US_MoteIn_Past_7_39_1140_RNASeq_R2.fq.gz
	
	$ sbatch KparkerRenamer_R_2.sh
	Submitted batch job 9693777

## 2022-03-28: Running TrimGalore

Create a list of all of the updated sample names.  

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$ls *R1.fq.gz > filenames.txt

Copy the filenames.txt list locally

	$ pwd 
	/c/Users/kpark/OneDrive/Documents/Barshis_Lab/2021-June_Mote/data/Sequencing/Porites_astreoides_fastq

	$  scp kpark049@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/filenames.txt ./
	kpark049@turing.hpc.odu.edu's password:
	filenames.txt                                                                                      100%   13KB 441.9KB/s   00:00

	$ head filenames.txt
	20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_34_180_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_34_360_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_34_420_RNASeq_R1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_34_750_RNASeq_R1.fq.gz

Open the filenames.txt in Notepad ++ and use regex to change it to the required format to run TrimGalore

	*Locally In Notepad ++*
	
	Find: (^\d.+\_)(\w)(\d)(.+)(\r)
	Replace: crun trim_galore --fastqc --paired \1\2\3\4 \1\22\4

File now looks like this:  

	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R2.fq.gz
	...

Use copy paste to separate it into 3 files with 73 samples (219 total samples/3 separate scripts) then add sbatch script header and additional lines to enable TrimGalore. 

	$ pwd  
	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/

	$ nano 2022-03-28_TrimGalore01.sh
	$ nano 2022-03-28_TrimGalore02.sh
	$ nano 2022-03-28_TrimGalore03.sh

Each 2022-03-28_TrimGalore0*.sh file looks like this with it's designated 73 lines of sample names:

	$ cat 2022-03-28_TrimGalore01.sh

	#!/bin/bash -l

	#SBATCH -o 2022-03-28_TrimGalore01.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=TrimGalore01

	enable_lmod

	module load container_env trim_galore

	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R2.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2.fq.gz
	...

	$ ls
	2022-03-28_TrimGalore01.sh  2022-03-28_TrimGalore02.sh  2022-03-28_TrimGalore03.sh  gunzip_test.py  KparkerRename.sh  renamer_advbioinf.py  renamer_KEP.py

Run each individual scrip to trim each sample's forward and reverse reads. 

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$ salloc
	salloc: Pending job allocation 9734126
	salloc: job 9734126 queued and waiting for resources
	salloc: job 9734126 has been allocated resources
	salloc: Granted job allocation 9734126
	salloc: Waiting for resource configuration
	salloc: Nodes coreV2-25-054 are ready for job
	This session will be terminated in 7 days. If your application requires a longer excution time, please use command "salloc -t N-0" where N is the number of days that you need.

	$ sbatch /cm/shared/courses/dbarshis/barshislab/KatieP/scripts/2022-03-28_TrimGalore01.sh
	Submitted batch job	 9734127

	$ sbatch /cm/shared/courses/dbarshis/barshislab/KatieP/scripts/2022-03-28_TrimGalore02.sh
	
	$ sbatch /cm/shared/courses/dbarshis/barshislab/KatieP/scripts/2022-03-28_TrimGalore03.sh

	$ sbatch /cm/shared/courses/dbarshis/barshislab/KatieP/scripts/2022-03-28_TrimGalore02.sh
	Submitted batch job 9734128

	$ sbatch /cm/shared/courses/dbarshis/barshislab/KatieP/scripts/2022-03-28_TrimGalore03.sh
	Submitted batch job 9734129

	$ squeue -u kpark049
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           9734129      main TrimGalo kpark049  R       1:14      1 coreV4-21-001
           9734128      main TrimGalo kpark049  R       1:24      1 coreV2-25-054
           9734127      main TrimGalo kpark049  R       3:57      1 coreV2-25-054
           9734126      main interact kpark049  R       7:08      1 coreV2-25-054
	
Jobs finished 2022-03-31, took a little under 3 days to complete.


## Notes on TrimGalore! 

--fastqc: Run FastQC in the default mode on the FastQ file once trimming is complete. 

--paired: performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files  
* To pass the validation test, both sequences of a sequence pair are required to have a certain minimum length which is governed by the option --length. If only one read passes this length threshold the other read can be rescued
* Using this option lets you discard too short read pairs without disturbing the sequence-by-sequence order of FastQ files which is required by many aligners. Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g. file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz

## 2022-04-12 Formatting Kenkel Porites transriptomes for STAR 

Moved trim reports into it's own directory

	$ /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/

	$ mkdir trim_reports 

	$ mv *_trimming_report.txt ./trim_reports

Create a new directory to move any Acropora cervicornis files. 

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons

	$ mkdir Acropora_cervicornis
	
	$ ls
	Acropora_cervicornis  Porites_astreoides

Move any A. cervicornis files in raw\_data\_fasqs to the new directory.

	**CAN'T FIND ACER FILES IN PAST DIRECTORY? SOMEWHERE BETWEEN RENAMING AND TRIMMING THEY DISSAPEARED??**

Create a "mapping" directory in raw\_data\_fasqs and "kenkelPast" directory 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/

	$ mkdir mapping 

	$ cd mapping/

	$ mkdir kenkelPast

	$ cd kenkelPast

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast


Copy Kenkel Genome to new directory and run MakeGenome script for Kenkel Past genome

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast
	
	$ scp /cm/shared/courses/dbarshis/barshislab/danb/taxons/Porites_astreoides/2021-12_MotePilotV1/Fastqs/mapping/kenkelPast/Porites_astreoides_LongestIsoform_suffixed.fasta ./

	$ ls  
	Porites_astreoides_LongestIsoform_suffixed.fasta

	$ nano MakeGenomeKenk.sh
	
	#!/bin/bash -l

	#SBATCH -o 2022-02-07_KenkelgenomeGenerate.txt
	#SBATCH -n 4
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=Kenkgenomegenerate

	enable_lmod

	module load container_env star

	STAR --runMode genomeGenerate --runThreadN 4 --genomeDir ./ --genomeFastaFiles Porites_astreoides_LongestIsoform_suffixed.fasta --genomeChrBinNbits 16

	$ salloc 
	salloc: Pending job allocation 9750125
	salloc: job 9750125 queued and waiting for resources
	salloc: job 9750125 has been allocated resources
	salloc: Granted job allocation 9750125

	$ sbatch MakeGenomeKenk.sh
	Submitted batch job 9750126

	$ squeue -u kpark049
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           9750126      main Kenkgeno kpark049  R       0:05      1 coreV1-22-012
           9750125      main interact kpark049  R       0:14      1 coreV1-22-018
	
After job finished (took ~ 3 hours)

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast
	$ ls
	2022-02-07_KenkelgenomeGenerate.txt  chrStart.txt          MakeGenomeKenk.sh                                      SAindex
	chrLength.txt                        Genome                Porites_astreoides_LongestIsoform_suffixed.fasta
	chrNameLength.txt                    genomeParameters.txt  Porites_astreoides_LongestIsoform_suffixed.fasta.save
	chrName.txt                          Log.out               SA

## 2022-04-14 Running MultiQC 

Run MultiQC 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$ salloc 
	salloc: Pending job allocation 9750557
	salloc: job 9750557 queued and waiting for resources
	salloc: job 9750557 has been allocated resources
	salloc: Granted job allocation 9750557
	
	$ module load container_env multiqc
	$ crun multiqc ./
	[WARNING]         multiqc : MultiQC Version v1.12 now available!
	[INFO   ]         multiqc : This is MultiQC v1.9
	[INFO   ]         multiqc : Template    : default
	[INFO   ]         multiqc : Searching   : /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs
	Searching 2223 files..  [####################################]  100%
	[INFO   ]        cutadapt : Found 438 reports
	[INFO   ]          fastqc : Found 438 reports
	[INFO   ]         multiqc : Compressing plot data
	[INFO   ]         multiqc : Report      : multiqc_report.html
	[INFO   ]         multiqc : Data        : multiqc_data
	[INFO   ]         multiqc : MultiQC complete

## 2022-04-19 Running STAR 
Run STAR scripts per genotype, testing on genotype 10 to start 

**THIS DID NOT WORK**  

I forgot to add in an ----outFileNamePrefix so all of the genotype 10 data was saved into one log file. All output from it was removed. 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$ nano STARMapKenkel_10.sh

	#!/bin/bash -l

	#SBATCH -o 2022-04-14_STARMapKenkel_10.txt
	#SBATCH -n 4
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=STARKenkel_10
	
	enable_lmod
	
	module load container_env star
	
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_180_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_360_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_420_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_750_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_180_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_360_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_750_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_1140_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_180_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_180_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_360_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_360_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_39_750_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_750_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_only_seq_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_only_seq_RNASeq_R2_val_2.fq.gz
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R2_val_2.fq.gz

	$ sbatch STARMapKenkel_10.sh
	Submitted batch job 9751644

## Notes on STAR
**--genomeDir:** specifies path to the directory (henceforth called ”genome directory” where the genome indices are stored

**--runThreadN 16:** option defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node

**--outSAMattributes:** a string of desired SAM attributes, in the order desired for the output SAM. Tags can be listed in any combination/order  
* **All** (argument includes the following)  
	* NH: number of loci the reads maps to: = 1 for unique mappers, > 1 for multimappers. Standard SAM tag.  
	* HI: multiple alignment index, starts with –outSAMattrIHstart (= 1 by default). Standard SAM tag.
	* AS: local alignment score, +1/ − 1 for matches/mismateches, score* penalties for indels and gaps. For PE reads, total score for two mates. Standard SAM tag.
	* NM: edit distance to the reference (number of mismatched + inserted + deleted bases) for each mate. Standard SAM tag.
	* nM: number of mismatches per (paired) alignment, not to be confused with NM, which is the number of mismatches+indels in each mate.
	* jM:B:c,M1,M2,... : intron motifs for all junctions (i.e. N in CIGAR): 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT. If splice junctions database is used, and a junction is annotated, 20 is added to its motif value.
	* MD: string encoding mismatched and deleted reference bases (see standard SAM specifications).Standard SAM tag.
	* jI:B:I,Start1,End1,Start2,End2,... : Start and End of introns for all junctions (1-based).
	* jM jI : attributes require samtools 0.1.18 or later, and were reported to be incompatible with some downstream tools such as Cufflinks.

**--genomeLoad:** mode of shared memory usage for the genome files. Only used with –runMode alignReads  
* **LoadAndRemove** = load genome into shared but remove it after run 

**--outFilterType:** type of filtering 
* **Normal** = standard filtering using only current alignment
 
**--outFilterMismatchNoverLmax:** alignment will be output only if its ratio of mismatches to *mapped*length is less than or equal to this value
* default 0.3

**--outSAMstrandField:** for unstranded RNA-seq data, generates required spliced alignments with XS strand attribute for Cufflinks/Cuffdiff
* XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed
* **intronMotif** = This option changes the output alignments: reads with inconsistent and/or non-canonical introns are filtered out

**--outFilterIntronMotifs:** filter alignment using their motifs, recommended to remove non-canonical junctions for Cufflinks runs    
* **RemoveNoncanonical** = filter out alignments that contain non-canonical junctions

**--outSAMtype BAM Unsorted:** output unsorted BAM file as Aligned.out.bam file  
* The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting. The order of the reads will match that of the input FASTQ(A)files only if one thread is used --runThread 1, and --outFilterType --BySJout is not used.
 
**--limitBAMsortRAM 5784458574:** int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option.

**--readFilesCommand:** command line to execute for each of the input file  
* **zcat** = to uncompress .gz files 

**--outReadsUnmapped:** output of unmapped and partially mapped (i.e. mapped only one mate of a paired end read) reads in separate file(s)
* **Fastx** = output in separate fasta/fastq files, Unmapped.out.mate1/2

**--outFilterMatchNminOverLread:** sam as outFilterMatchNmin, but normalized to the read length (sum of mates’ lengths for paired-end reads)  
* default: 0.66

**--outFilterScoreMinOverLread:** alignment will be output only if its score is higher than or equal to this value but normalized to read length (sum of mates’ lengths for paired-end reads)  
* default 0.66

**--readFilesIn:** name(s) (with path) of the files containing the sequences to be mapped (e.g. RNA-seq FASTQ files)  
* For paired-end reads, use comma separated list for read1, followed by space, followed by comma separated list for read2

**--outFileNamePrefix:** string, output files name prefix (including full or relative path). Can only be defined on the command line

## 2022-04-20 Retesting running STAR on one sample  
Retrying STAR using one sample "20210608T1400\_CBASS\_US\_MoteIn\_Past\_10\_30\_1140" with the outputfile command 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastq  

	$ nano STARMapKenkel_test.sh
	#!/bin/bash -l

	#SBATCH -o 2022-04-20_STARMapKenkel_test.txt
	#SBATCH -n 4
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=STARKenkel_test
	
	enable_lmod
	
	module load container_env star
	
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1_val_1.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix ${i%_R1_val_1.fq.gz}_2kenkel

	$ salloc 
	salloc: Pending job allocation 9752138
	salloc: job 9752138 queued and waiting for resources
	salloc: job 9752138 has been allocated resources
	salloc: Granted job allocation 9752138
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.  

	$ sbatch STARMapKenkel_test.sh
	Submitted batch job 9752139

	$ squeue -u kpark049
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           9752139      main STARKenk kpark049  R       0:12      1 coreV4-21-001
           9752138      main interact kpark049  R       0:26      1 coreV1-22-001

This worked as needed, but need to put R1 and R2 on the same line to have the mapping in one file output. 

## 2022-04-21 Finally Running STAR!

Started off with list of sample names from TrimGalore

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fast

	$ ls *R1_val_1.fq.gz > trimgalore_filenames.txt

	$ head trimgalore_filnames.txt
	20210608T1400_CBASS_US_MoteIn_Past_10_30_1140_RNASeq_R1_val_1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_180_RNASeq_R1_val_1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1_val_1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R1_val_1.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_750_RNASeq_R1_val_1.fq.gz

Safe copied trimgalore_filenames.txt to local and edited in Notpad ++

	find: ^(\d+.+RNASeq\_)(R)(1)(_val_)(1)(.fq.gz)
	replace: STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn \1\2\3\4\5\6 \1\22\42\6 --outFileNamePrefix \12kenkel_

Then moved genotypes in groups of two to separate files to make scripts and copied to cluster

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fast

	$ nano STARMapKenkel_01_02.sh

	#!/bin/bash -l

	#SBATCH -o 2022-04-21_STARMapKenkel_01_02.txt
	#SBATCH -n 4
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=STARKenkel_01_02
	
	enable_lmod
	
	module load container_env star
	
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_30_1140_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_30_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_30_1140_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_30_180_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_30_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_30_180_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_30_360_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_30_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_30_360_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_30_420_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_30_420_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_30_420_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_30_750_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_30_750_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_30_750_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_34_1140_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_34_1140_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_34_1140_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_34_180_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_34_180_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_34_180_RNASeq_2kenkel_
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210610T1400_CBASS_US_MoteOff_Past_1_34_360_RNASeq_R1_val_1.fq.gz 20210610T1400_CBASS_US_MoteOff_Past_1_34_360_RNASeq_R2_val_2.fq.gz --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_1_34_360_RNASeq_2kenkel_
	...

Did the same for geonotype 3 and 4, 5 and 6, 7 and 8, and 9 and 10. Ran all on turing. 

	$ salloc 
	salloc: Pending job allocation 9753159
	salloc: job 9753159 queued and waiting for resources
	salloc: job 9753159 has been allocated resources
	salloc: Granted job allocation 9753159
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.

	$ sbatch STARMapKenkel_01_02.sh
	Submitted batch job 9753160  

	$ sbatch STARMapKenel_03_04.sh
	Submitted batch job 9753162

	$ sbatch STARMapKenkel_05_06.sh
	Submitted batch job 9753163  

	$ sbatch STARMapKenkel_07_08.sh
	Submitted batch job 9753166

	$ sbatch STARMapKenkel_09_10.sh
	Submitted batch job 9753165

	$ squeue -u kpark049
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           9753166      main STARKenk kpark049  R       0:10      1 coreV2-25-007
           9753165      main STARKenk kpark049  R       1:53      1 coreV1-22-012
           9753163      main STARKenk kpark049  R       2:14      1 coreV2-22-024
           9753162      main STARKenk kpark049  R       2:47      1 coreV2-22-010
           9753160      main STARKenk kpark049  R       3:42      1 coreV2-25-knc-010
           9753161      main interact kpark049  R       3:17      1 coreV1-22-018
           9753159      main interact kpark049  R       4:39      1 coreV1-22-018
	
	
Started jobs at 3:05 PM
Last job ended at 7:30 AM next day (04/22/2022)

	Slurm Job_id=9753166 Name=STARKenkel_07_08 Ended, Run time 14:02:38, COMPLETED, ExitCode 0  
	Slurm Job_id=9753160 Name=STARKenkel_01_02 Ended, Run time 14:25:19, COMPLETED, ExitCode 0   
	3 of 1,578
	Slurm Job_id=9753162 Name=STARKenkel_03_04 Ended, Run time 14:33:21, COMPLETED, ExitCode 0 
	Slurm Job_id=9753163 Name=STARKenkel_05_06 Ended, Run time 15:49:59, COMPLETED, ExitCode 0
	Slurm Job_id=9753165 Name=STARKenkel_09_10 Ended, Run time 16:23:17, COMPLETED, ExitCode 0

## Running multiqc on trimgalore and STAR output

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs

	$ salloc 
	salloc: Pending job allocation 9754717
	salloc: job 9754717 queued and waiting for resources
	salloc: job 9754717 has been allocated resources
	salloc: Granted job allocation 9754717  

	$ enable_lmod	
	
	$ module load container_env multiqc
	
	$ crun multiqc ./

	[WARNING]         multiqc : MultiQC Version v1.12 now available!
	[INFO   ]         multiqc : This is MultiQC v1.9
	[INFO   ]         multiqc : Template    : default
	[INFO   ]         multiqc : Searching   : /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs
	Searching 3908 files..  [####################################]  100%
	[INFO   ]            star : Found 215 reports
	[INFO   ]        cutadapt : Found 438 reports
	[INFO   ]          fastqc : Found 438 reports
	[INFO   ]         multiqc : Compressing plot data
	[WARNING]         multiqc : Previous MultiQC output found! Adjusting filenames..
	[WARNING]         multiqc : Use -f or --force to overwrite existing reports instead
	[INFO   ]         multiqc : Report      : multiqc_report_1.html
	[INFO   ]         multiqc : Data        : multiqc_data_1
	[INFO   ]         multiqc : MultiQC complete  

## 2022-05-06 Copying re-sequenced files into a new directory

Making new directory to analyze re-sequenced files  
  
	$ pwd 
	/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/

	$ mkdir reseq_raw_data_fastqs 

	$ ls 
	KPunzip.txt  old              reseq_raw_data_fastqs  X204SC21081158-Z01-F002
	__MACOSX     raw_data_fastqs  unzip_raw.sh           X204SC21081158-Z01-F002.zip

Copying files to new directory 

	* immediately after cluster login * 
	$ cp /RC/group/rc_barshis_lab/taxonarchive/Porites_astreoides/2022-04-27_MotePilotResequencing/*.zip /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs 

Check file copy 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs

	$ ls 
	X204SC21081158-Z01-F012_01.zip	X204SC21081158-Z01-F012_02.zip

## 2022-05-08 Running md5cheksum  

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs  

	$ salloc 
	salloc: Pending job allocation 9756952
	salloc: job 9756952 queued and waiting for resources
	salloc: job 9756952 has been allocated resources
	salloc: Granted job allocation 9756952

	
	$ nano ChekSum.sh 
	
	#!/bin/bash -l

	#SBATCH -o 2022-05-08_md5.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=md5
	
	md5sum X204SC21081158-Z01-F012_01.zip  
	md5sum X204SC21081158-Z01-F012_02.zip  

	$ sbatch ChekSum.sh
	Submitted batch job 9756953

	$ squeue -u kpark049
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           9756953      main      md5 kpark049  R       0:09      1 coreV2-22-005
           9756952      main interact kpark049  R       0:41      1 coreV2-22-005 

	$ cat 2022-05-08_md5.txt
	12ab39a924ded5216df7deabcf255a79  X204SC21081158-Z01-F012_01.zip
	7a45743a0891fcaa2c4e383439e03308  X204SC21081158-Z01-F012_02.zip

Running unzip script 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs
	
	$ nano reseq_unzip_raw.sh
	#!/bin/bash -l

	#SBATCH -o KPunzip.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=KPunzip
	
	unzip X204SC21081158-Z01-F012_01.zip  
	unzip X204SC21081158-Z01-F012_02.zip
	
	$ salloc
	salloc: Pending job allocation 9756958
	salloc: job 9756958 queued and waiting for resources
	salloc: job 9756958 has been allocated resources
	salloc: Granted job allocation 9756958
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.

	$ sbatch reseq_unzip_raw.sh
	Submitted batch job 9756959

	$ squeue -u kpark049
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           9756959      main reseq_un kpark049  R       0:07      1 coreV2-22-005
           9756958      main interact kpark049  R       0:39      1 coreV2-22-005  

## 2022-05-09: Copying unzipped files into main directory, running renamer 

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs/X204SC21081158-Z01-F012_01/raw_data

	
	$ ls
	R456  R46   R48  R52  R56  R61  R66  R70  R74  R8   R83  R91                 RE152  RE157  RE186  RE192  RR147  RR158
	R457  R460  R49  R53  R57  R63  R67  R71  R75  R80  R85  R96                 RE153  RE158  RE187  RE195  RR153  RR161
	R458  R461  R50  R54  R58  R64  R69  R72  R76  R81  R89  Rawdata_Readme.pdf  RE154  RE183  RE188  RE200  RR155  RR163
	R459  R47   R51  R55  R59  R65  R7   R73  R79  R82  R90  RE151               RE156  RE185  RE189  RR143  RR157  RR164


	$ cp ./R*/*.fq.gz /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs/

Did the same thing for X204SC21081158-Z01-F012_02

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs

	$ ls 
	2022-05-08_md5.txt  R283r_1.fq.gz  R318_1.fq.gz   R35_1.fq.gz    R4_1.fq.gz    R57_1.fq.gz  RE151_1.fq.gz
	ChekSum.sh          R283r_2.fq.gz  R318_2.fq.gz   R351r_1.fq.gz  R42_1.fq.gz   R57_2.fq.gz  RE151_2.fq.gz
	R100_1.fq.gz        R284_1.fq.gz   R319_1.fq.gz   R351r_2.fq.gz  R42_2.fq.gz   R58_1.fq.gz  RE152_1.fq.gz
	R100_2.fq.gz        R284_2.fq.gz   R319_2.fq.gz   R352_1.fq.gz   R423_1.fq.gz  R58_2.fq.gz  RE152_2.fq.gz
	R101_1.fq.gz        R285_1.fq.gz   R3_1.fq.gz     R352_2.fq.gz   R423_2.fq.gz  R59_1.fq.gz  RE153_1.fq.gz
	R101_2.fq.gz        R285_2.fq.gz   R320_1.fq.gz   R35_2.fq.gz    R429_1.fq.gz  R59_2.fq.gz  RE153_2.fq.gz
	R102_1.fq.gz        R286r_1.fq.gz  R320_2.fq.gz   R355_1.fq.gz   R429_2.fq.gz  R61_1.fq.gz  RE154_1.fq.gz
	R102_2.fq.gz        R286r_2.fq.gz  R32_1.fq.gz    R355_2.fq.gz   R4_2.fq.gz    R61_2.fq.gz  RE154_2.fq.gz
	R103_1.fq.gz        R287_1.fq.gz   R32_2.fq.gz    R356_1.fq.gz   R43_1.fq.gz   R63_1.fq.gz  RE156_1.fq.gz
	R103_2.fq.gz        R287_2.fq.gz   R322r_1.fq.gz  R356_2.fq.gz   R43_2.fq.gz   R63_2.fq.gz  RE156_2.fq.gz
	R104_1.fq.gz        R288_1.fq.gz   R322r_2.fq.gz  R36_1.fq.gz    R437_1.fq.gz  R64_1.fq.gz  RE157_1.fq.gz
	R104_2.fq.gz        R288_2.fq.gz   R323_1.fq.gz   R36_2.fq.gz    R437_2.fq.gz  R64_2.fq.gz  RE157_2.fq.gz
	R105_1.fq.gz        R289_1.fq.gz   R323_2.fq.gz   R364_1.fq.gz   R438_1.fq.gz  R65_1.fq.gz  RE158_1.fq.gz
	R105_2.fq.gz        R289_2.fq.gz   R324_1.fq.gz   R364_2.fq.gz   R438_2.fq.gz  R65_2.fq.gz  RE158_2.fq.gz
	...

Running Renamer: used original "R\_1\_fastq\_rename.txt" file and used find and replace to add in "reseq" to new file name. 

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs

	$ nano R_1_fastq_reseq_rename.txt
	OldName	NewName
	R100_1.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1_reseq.fq.gz
	R101_1.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R1_reseq.fq.gz
	R10_1.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R1_reseq.fq.gz
	R102_1.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R1_reseq.fq.gz
	R103_1.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R1_reseq.fq.gz
	...

	$ nano R_2_fastq_reseq_rename.txt
	OldName	NewName
	R100_2.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2_reseq.fq.gz
	R101_2.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R2_reseq.fq.gz
	R10_2.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R2_reseq.fq.gz
	R102_2.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R2_reseq.fq.gz
	R103_2.fq.gz	20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R2_reseq.fq.gz
	...

	$ nano  R_1_fastq_reseq_rename.sh
	#!/bin/bash -l

	#SBATCH -o reseq_fastq_rename_R1.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=reseqFastqR1
	
	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_1_fastq_reseq_rename.txt

	$ nano  R_2_fastq_reseq_rename.sh
	#!/bin/bash -l

	#SBATCH -o reseq_fastq_rename_R2.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=reseqFastqR2
	
	/cm/shared/courses/dbarshis/barshislab/KatieP/scripts/renamer_KEP.py R_2_fastq_reseq_rename.txt

	$ salloc 
	
	$ sbatch R_1_fastq_reseq_rename.sh
	
	$ sbatch R_2_fastq_reseq_rename.sh

	$ ls
	20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R1_reseq.fq.gz    20210610T1400_CBASS_US_MoteOff_Past_2_37_420_RNASeq_R2_reseq.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_360_RNASeq_R2_reseq.fq.gz    20210610T1400_CBASS_US_MoteOff_Past_2_39_1140_RNASeq_R1_reseq.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R1_reseq.fq.gz    20210610T1400_CBASS_US_MoteOff_Past_2_39_1140_RNASeq_R2_reseq.fq.gz
	20210608T1400_CBASS_US_MoteIn_Past_10_30_420_RNASeq_R2_reseq.fq.gz    20210610T1400_CBASS_US_MoteOff_Past_2_39_360_RNASeq_R1_reseq.fq.gz
	...

## 2022-05-19: Running TrimGalore on reseq samples  

	$ pwd /cm/shared/courses/dbarshis/barshislab/KatieP/scripts/

	$ nano 2022-05-19_trimgalore_reseq01.sh
	#!/bin/bash -l

	#SBATCH -o 2022-05-19_trimgalore_reseq01.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=TrimGalore01
	
	enable_lmod
	
	module load container_env trim_galore
	
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_6_39_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_6_39_1140_RNASeq_R2_reseq.fq.gz
	...

	$ nano 2022-05-19_trimgalore_reseq02.sh
	#!/bin/bash -l
	
	#SBATCH -o 2022-05-19_trimgalore_reseq02.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=TrimGalore02
	
	enable_lmod
	
	module load container_env trim_galore
	
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_9_37_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_6_39_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_6_39_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_7_39_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_7_39_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_8_39_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_8_39_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_9_39_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_39_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_39_420_RNASeq_R2_reseq.fq.gz
	crun trim_galore --fastqc --paired 20210608T1400_CBASS_US_MoteIn_Past_6_30_750_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_6_30_750_RNASeq_R2_reseq.fq.gz
	...

Running reseq TrimGalore scripts 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs  
	
	$ salloc

	$ sbatch 2022-05-19_trimgalore_reseq01.sh

	$ sbatch 2022-05-19_trimgalore_reseq02.sh


## 2022-05-23: MultiQC and STAR 

Running MultiQC

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs

	$ $ salloc
	salloc: Pending job allocation 9761166
	salloc: job 9761166 queued and waiting for resources
	salloc: job 9761166 has been allocated resources
	salloc: Granted job allocation 9761166
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.
	
	$ module load container_env multiqc

	$ crun multiqc ./

	$ [WARNING]         multiqc : MultiQC Version v1.12 now available!
	[INFO   ]         multiqc : This is MultiQC v1.9
	[INFO   ]         multiqc : Template    : default
	[INFO   ]         multiqc : Searching   : /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs
	Searching 2178 files..  [####################################]  100%
	[INFO   ]        cutadapt : Found 312 reports
	[INFO   ]          fastqc : Found 312 reports
	[INFO   ]         multiqc : Compressing plot data
	[INFO   ]         multiqc : Report      : multiqc_report.html
	[INFO   ]         multiqc : Data        : multiqc_data
	[INFO   ]         multiqc : MultiQC complete

Running STAR -upgraded 22-05-19_trimgalore_reseq0*.sh files with text editor to STAR format, using same genome location from previous STAR run 

	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs

	$ salloc
	salloc: Pending job allocation 9761167
	salloc: job 9761167 queued and waiting for resources
	salloc: job 9761167 has been allocated resources
	salloc: Granted job allocation 9761167
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.  

	$ nano STARKenkel_reseq01.sh
	#!/bin/bash -l

	#SBATCH -o 2022-05-23-STARKenkel_reseq01.txt
	#SBATCH -n 4
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=STARKenkel_reseq01
	
	enable_lmod
	
	module load container_env star
	
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_34_1140_RNASeq_reseq_2kenkel
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_6_37_1140_RNASeq_reseq_2kenkel
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_T0_0_RNASeq_reseq_2kenkel
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_7_37_1140_RNASeq_reseq_2kenkel
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_8_37_1140_RNASeq_reseq_2kenkel
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_9_37_1140_RNASeq_reseq_2kenkel
	STAR --genomeDir /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/ --runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2 --readFilesIn 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R1_reseq.fq.gz 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_R2_reseq.fq.gz --outFileNamePrefix 20210608T1400_CBASS_US_MoteIn_Past_10_37_1140_RNASeq_reseq_2kenkel
	...

Did the same for 2022-05-23-STARKenkel\_reseq02.txt

	$ sbatch STARKenkel_reseq01.sh
	Submitted batch job 9761168

	$ sbatch STARKenkel_reseq02.sh
	Submitted batch job 9761169

## 2022-05-25: Running MultiQC on TrimGalore and STAR output 
	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs
	
	$ salloc
	salloc: Pending job allocation 9761397
	salloc: job 9761397 queued and waiting for resources
	salloc: job 9761397 has been allocated resources
	salloc: Granted job allocation 9761397
	This session will be terminated in 7 days. If your application requires
	a longer excution time, please use command "salloc -t N-0" where N is the
	number of days that you need.
	
	$ enable_lmod
	
	$ module load container_env multiqc

	$ crun multiqc ./

	Searching 4955 files..  [####################################]  100%
	[INFO   ]            star : Found 396 reports
	[INFO   ]        cutadapt : Found 312 reports
	[INFO   ]          fastqc : Found 312 reports
	[INFO   ]         multiqc : Compressing plot data
	[WARNING]         multiqc : Previous MultiQC output found! Adjusting filenames..
	[WARNING]         multiqc : Use -f or --force to overwrite existing reports instead
	[INFO   ]         multiqc : Report      : multiqc_report_1.html
	[INFO   ]         multiqc : Data        : multiqc_data_1
	[INFO   ]         multiqc : MultiQC complete  

Copying output files locally. 

	$ pwd 
	/c/Users/kpark/OneDrive/Documents/Barshis_Lab/2021-June_Mote/data/Sequencing

	$ scp kpark049@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs/multiqc_report*.html ./
	kpark049@turing.hpc.odu.edu's password:
	multiqc_report.html                                        100% 6853KB   6.5MB/s   00:01
	multiqc_report_1.html                                      100% 8321KB   5.3MB/s   00:01  

Renamed locally as: 

multiqc\_report.html ---> 2022-05-23\_reseq\_multiqc\_report\_TG.html  
multiqc\_report_1.html  ----> 2022-05-25\_reseq\_multiqc\_report\_TG\_STAR.html

## 2022-07-12 Salmon Test 01

Move STAR BAM files from reseq'd files to directory with other BAM files   
	
	$ pwd  
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/reseq_raw_data_fastqs

	$ mv *_reseq_2kenkelAligned.out.bam /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/bam_files 
	
	$ pwd 
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/bam_files 
	
	$ ls
	20210610T1400_CBASS_US_MoteOff_Past_3_34_360_RNASeq_reseq_2kenkelAligned.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_3_34_420_RNASeq_2kenkel_Aligned.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_3_34_420_RNASeq_reseq_2kenkelAligned.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_3_34_750_RNASeq_2kenkel_Aligned.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_3_34_750_RNASeq_reseq_2kenkelAligned.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_3_37_1140_RNASeq_2kenkel_Aligned.out.bam
	...  

	$ nano 2022-07-12_SalmonTest01.sh

	#!/bin/bash -l

	#SBATCH -o 2022-07-12_SalmonTest01.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=SalmonTest01
	
	enable_lmod
	
	module load container_env salmon/1.9.0
	
	crun salmon quant -t /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/Porites_astreoides_LongestIsoform_suffixed.fasta <LIBTYPE> -a ./20210610T1400_CBASS_US_MoteOff_Past_3_30_750_RNASeq_2kenkel_Aligned.out.bam -o salmon_quant

This did not work, because I did not specify the LIBTYPE parameter. 

	$ head 2022-07-12_SalmonTest01.txt
	/home/kpark049/.bash_profile: line 5: /home/kpark049/.turing_bash_profile: No such file or directory 
	/var/spool/slurmd/job9781846/slurm_script: line 13: LIBTYPE: No such file or directory

##2022-07-14 Salmon Test 02

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/bam_files

	$ nano 2022-07-14_SalmonTest02.sh

	#!/bin/bash -l
	
	#SBATCH -o 2022-07-14_SalmonTest02.txt
	#SBATCH -n 1
	#SBATCH --mail-user=kpark049@odu.edu
	#SBATCH --mail-type=END
	#SBATCH --job-name=SalmonTest02
	
	enable_lmod
	
	module load container_env salmon/1.9.0
	
	crun salmon quant -t /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/Porites_astreoides_LongestIsoform_suffixed.fasta --libType A -a 20210610T1400_CBASS_US_MoteOff_Past_3_30_750_RNASeq_2kenkel_Aligned.out.bam -o salmon_quant


It worked!! 

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/bam_files/salmon_quant

	$ ls
	aux_info  cmd_info.json  libParams  logs  quant.sf

	$ head quant.sf
	Name    Length  EffectiveLength TPM     NumReads
	GCKDGN101CF7JK_Past     223     68.354  -nan    0.000
	GCKDGN101CAZ1A_Past     363     148.513 -nan    0.000
	GCKDGN101ANI4S_Past     214     70.510  -nan    0.000
	GCKDGN101BV6U3_Past     277     80.600  -nan    0.000
	GCKDGN101BQG8Q_Past     412     196.116 -nan    0.000
	GCKDGN101CG9W2_Past     255     71.884  -nan    0.980
	GCKDGN101A7H4U_Past     352     138.399 -nan    0.000
	GCKDGN101ATRIP_Past     304     97.932  -nan    0.000
	GCKDGN101CINM3_Past     337     125.002 -nan    0.000

2022-08 (Some info lost due to markdown save error)

Combined quant.sf files from salmon output for reseq files using following python script. Input was a .txt file list with the original and reseq file names separated by tabs

	#!/usr/bin/env python
	
	import sys
	#sys.argv[1] file with list of your paired files to combine
	
	columntoextract = 4
	
	def make_dict1(file):
		fin = open(file, 'r')
		dict={}
		count=0
		for line in fin:
			count+=1
			line=line.rstrip()
			cols=line.split('\t') #for tab-delimited text files
			if count > 1:
				dict[cols[0]]=float(cols[columntoextract])
		fin.close()
		print('Read in first file')
		return dict
		
	InfileList=open(sys.argv[1], 'r')
	FileCount=0
	for line in InfileList:
		FileCount+=1
		line=line.rstrip()
		Files=line.split('\t')
		FirstFileDict=make_dict1(Files[0])
		SecondFile=open(Files[1], 'r')
		SecondFileLine=0
		Outfile=open('%s_merged.txt'%(Files[0][:-3]), 'w')
		for Record in SecondFile:
			SecondFileLine+=1
			Record=Record.rstrip()
			Items=Record.split('\t')
			if SecondFileLine==1:
				Outfile.write('ContigName\tNumReads\r')
			if SecondFileLine>1:
				Newcount=float(Items[columntoextract])+FirstFileDict[Items[0]]
				Outfile.write('%s\t%.3f\r'%(Items[0],Newcount))
		SecondFile.close()
		Outfile.close()
		print('Merged %s and %s into %s' %(Files[0], Files[1],'%s_merged.txt'%(Files[0][:-3])))

Read counts for files separated into read counts for individual contigs. DEGs calculated with sig value of 0.10 

	#!/usr/bin/env python
	
	Usage = """
	DEcalc.py - version 1.0
	Calculate the number of Differentially Expressed Genes for each contrast 
	and the up or down regulation of DEGs.
	Usage:
		DEcalc.py *.txt
	"""
	
	# import modules 
	import sys # lets us access system-specific parameters and functions
	
	sig = float(0.10)
	
	if len(sys.argv)<1:
		print(Usage) #print usage statements to the screen
	else:
		FileList = sys.argv[1:]
		FileNumber = 0
		CombinedOutFileName = "2022-11-01_dedup_test.tab"
		CombinedOutFile = open(CombinedOutFileName, 'w')
		for InFileName in FileList:
			InFile = open(InFileName, 'r') # open the infile
			Contrast = InFileName[:-4]
			Headerline = "Contrast,NumDEGsUp,NumDEGsDown,TotalNumDEGs"
			if FileNumber == 0:
				CombinedOutFile.write(Headerline + '\n')
			LineNumber = 0
			TotalNumDEGs = 0
			NumDEGsUp = 0
			NumDEGsDown = 0
			for Line in InFile:
				LineNumber+=1 # Line Number + 1
				if LineNumber > 1 :
					Line = Line.strip("\n").strip("\r")
					List = Line.split('\t')
					Contig = List[0]
					baseMean = List[1]
					log2FoldChange = float(List[2])
					lfcSE = List[3]
					stat = List[4]
					if List[5] == 'NA' or List[6] == 'NA':
						continue 
					else:
						pvalue = float(List[5])
						padj = float(List[6])
					
					# "Contrast,NumDEGsUp,NumDEGsDown,TotalNumDEGs"
					if padj < sig:
						TotalNumDEGs += 1
						if log2FoldChange > 0:
							NumDEGsUp += 1
						elif log2FoldChange < 0:
							NumDEGsDown += 1
	
			CombinedOutFile.write('%s\t%d\t%d\t%d\n' %(Contrast,NumDEGsUp,NumDEGsDown,TotalNumDEGs))
			InFile.close()
			FileNumber+= 1
			print("finished processing file: %s"%(InFileName))

 
		 

2022-11-01
Testing de-duplication through STAR

	$ pwd
	/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/bam_files  

	$ mkdir Star_bam_dedup_sandbox

	$ cp ./20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned.out.bam ./Star_bam_dedup_sandbox
	 
	$ cd Star_bam_dedup_sandbox

	$ ls 
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned.out.bam

	$ enable_lmod

	$ module load container_env star

	$ STAR --runThreadN 32 --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType UniqueIdenticalNotMulti --inputBAMfile 20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned.out.bam --outFileNamePrefix 20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned_dedup

	$ ls
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned_dedupLog.out
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned_dedupProcessed.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned.out.bam

	$ enable_lmod

	$ module load container_env salmon/1.9.0

	crun salmon quant -t /cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/mapping/kenkelPast/Porites_astreoides_LongestIsoform_suffixed.fasta --libType A -a 20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned_dedupProcessed.out.bam -o salmon_quant 

	$ ls
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned_dedupLog.out
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned_dedupProcessed.out.bam
	20210610T1400_CBASS_US_MoteOff_Past_5_T0_0_RNASeq_2kenkel_Aligned.out.bam
	salmon_quant/

	$ cd salmon_quant/
	aux_info  cmd_info.json libParams  logs  quant.sf

	$ cp ./quant.sf ./dedup_test_quant.sf


On local machine:  

	$ scp kpark049@turing.hpc.odu.edu:/cm/shared/courses/dbarshis/barshislab/KatieP/taxons/Porites_astreoides/2021-12_MotePilotV1/raw_data_fastqs/bam_files/STAR_bam_dedup_sandbox/salmon_quant/dedup_test_quant.sf ./

	$ head dedup_test_quant.sf

	$  Name    Length  EffectiveLength TPM     NumReads
	GCKDGN101CF7JK_Past     223     36.362  -nan    0.000
	GCKDGN101CAZ1A_Past     363     133.229 -nan    0.000
	GCKDGN101ANI4S_Past     214     34.489  -nan    70.000
	GCKDGN101BV6U3_Past     277     61.958  -nan    13.000
	GCKDGN101BQG8Q_Past     412     180.996 -nan    2.000
	GCKDGN101CG9W2_Past     255     49.726  -nan    14.000
	GCKDGN101A7H4U_Past     352     122.971 -nan    20.000
	GCKDGN101ATRIP_Past     304     81.476  -nan    10.000
	GCKDGN101CINM3_Past     337     109.455 -nan    1.000  



	


	
	

