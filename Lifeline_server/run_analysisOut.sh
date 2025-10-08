#!/bin/bash
#SBATCH --job-name=MAmestoy
#SBATCH --time=24:11:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=20000
#SBATCH --export=NONE
#SBATCH --get-user-env=L60

pwd
module load RPlus/4.2.1-foss-2022a-v22.12.1
#1500 + 1:30
#Rscript management_extraction.r
Rscript dev_analyse_out.R MainBMI
