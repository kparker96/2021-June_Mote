#!/bin/bash -l

#SBATCH -o 2022-05-08_md5.txt
#SBATCH -n 1
#SBATCH --mail-user=kpark049@odu.edu
#SBATCH --mail-type=END
#SBATCH --job-name=md5

md5sum X204SC21081158-Z01-F012_01.zip  
md5sum X204SC21081158-Z01-F012_02.zip
