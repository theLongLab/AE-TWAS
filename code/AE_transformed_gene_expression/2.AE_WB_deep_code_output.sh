#!/bin/bash
#SBATCH --job-name=34_D
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-0:0:0
#SBATCH --partition=theia,mtst,cpu2022,cpu2023

#python AE_WB_deep.py $1
python 2.AE_WB_deep_code_output.py "34" "/work/long_lab/qli/WGCNA/GTEX/12K_TPM1_1-22/Trial3a_20201220/34.L7.pth"
#python 2.AE_WB_deep_code_output.py 72 /work/long_lab/qli/WGCNA/GTEX/12K_TPM1_1-22/Trial3a_20201220/72.L7.pth
