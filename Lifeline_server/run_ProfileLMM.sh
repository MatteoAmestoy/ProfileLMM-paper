#!/bin/bash
#SBATCH --job-name=MAmestoy
#SBATCH --time=36:11:30
#SBATCH --cpus-per-task=4
#SBATCH --mem=20000
#SBATCH --export=NONE
#SBATCH --get-user-env=L60

pwd
module load RPlus/4.2.1-foss-2022a-v22.12.1
#1500 + 1:30
#Rscript management_extraction.r
Rscript main_ProfileLMM.R 2 60000 60 10000 5000 'DiasStab'
Rscript dev_analyse_out.R DiasStab 10000 8000
Rscript dev_clusterPostprocess.R DiasStab
