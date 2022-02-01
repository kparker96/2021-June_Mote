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


