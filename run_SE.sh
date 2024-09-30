#!/bin/bash 

#SBATCH --job-name=SE_run
#SBATCH --output=/home/ndo/slurm_script/SE_run%j.out
#SBATCH --error=/home/ndo/slurm_script/SE_run%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=4G

eval "$(conda shell.bash hook)"
eval "$(/home/ndo/miniconda3/bin/conda shell.bash hook)"
conda activate /home/ndo/miniconda3/envs/nextflow

cd ~/nextflow_SE/
nextflow run SE.nf -resume
