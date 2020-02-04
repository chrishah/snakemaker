#!/bin/bash
#
#SBATCH -J sm
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --mem=50G
#SBATCH --time=10:00:00

#SBATCH --qos normal_0064
#SBATCH --partition=mem_0064

#SBATCH --output slurm-%j.out
#SBATCH --error slurm-%j.err

#From here on out the magic happens
echo -e "\n[ $(date) ]\n"
echo "Running under shell '${SHELL}' in directory '$(pwd)' using $SLURM_NTASKS_PER_NODE tasks per node on $SLURM_JOB_NODES nodes ($SLURM_NODEID)"
echo "Host: $HOSTNAME"
echo "Job: $SLURM_JOB_NAME"

#load singularity
module load go/1.11 singularity/3.4.1

#activate snakemake environment
source ~/.bashrc
conda activate snakemake


echo -e "\n$(date) - Started\n"

#run snakemake
snakemake -p -j 16 --use-singularity

echo -e "\n$(date) - Finished\n"
