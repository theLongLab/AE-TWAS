#!/bin/bash
#SBATCH --job-name=AETWAS_S
#SBATCH --error=%x-%j.error
#SBATCH --out=%x-%j.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-0:0:0
#SBATCH --partition=gpu-v100
#SBATCH --gres=gpu:1

#python AE_WB_deep.py $1
python AE_WB_shallow.py $1
