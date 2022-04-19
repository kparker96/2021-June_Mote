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


Run STAR scripts per genotype, testing on genotype 10 to start 

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