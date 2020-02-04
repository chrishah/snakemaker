#!/bin/bash
#
#SBATCH -J sm
#SBATCH -N 1
#SBATCH --ntasks-per-node=30
#SBATCH --ntasks-per-core=1
#SBATCH --mem=50G
#SBATCH --time=07:00:00

#SBATCH --qos mem_0096
#SBATCH --partition=mem_0096

#SBATCH --output slurm-%j.out
#SBATCH --error slurm-%j.err

#From here on out the magic happens
echo -e "\n[ $(date) ]\n"
echo "Running under shell '${SHELL}' in directory '$(pwd)' using $SLURM_NTASKS_PER_NODE tasks per node on $SLURM_JOB_NODES nodes ($SLURM_NODEID)"
echo "Host: $HOSTNAME"
echo "Job: $SLURM_JOB_NAME"

#activate snakemake environment
source ~/.bashrc
conda activate snakemake


echo -e "\n$(date) - Started\n"

#run snakemake
snakemake -p --use-singularity

echo -e "\n$(date) - Finished\n"
