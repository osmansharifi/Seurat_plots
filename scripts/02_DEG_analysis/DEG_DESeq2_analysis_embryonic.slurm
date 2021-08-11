#!/bin/bash
#
#SBATCH --mail-user=vhaghani@ucdavis.edu                                                # User email to receive updates
#SBATCH --mail-type=ALL                                                                 # Get an email when the job begins, ends, or if it fails
#SBATCH -p production                                                                   # Partition, or queue, to assign to
#SBATCH -J embryonic_deseq2                                                             # Name for job
#SBATCH -o embryonic_deseq2_slurm.j%j.out                                               # File to write STDOUT to
#SBATCH -e embryonic_deseq2_slurm.j%j.err                                               # File to write error output to
#SBATCH -N 1                                                                            # Number of nodes/computers
#SBATCH -n 10                                                                           # Number of cores
#SBATCH -t 7:00:00                                                                      # Ask for no more than 7 hours
#SBATCH --mem=16gb                                                                      # Ask for no more than 16 GB of memory
#SBATCH --chdir=/share/korflab/home/viki/snRNA-seq-pipeline/scripts/02_DEG_analysis     # Directory I want the job to run in

# Run aklog to deal with SLURM bug
aklog

# Source .profile
source .profile

# Load R
module load R

# Fail on weird errors
set -o nounset
set -o errexit
set -x

# Run the snakemake!
Rscript DEG_DESeq2_analysis_embryonic.R

# Print out various information about the job
env | grep SLURM                                               # Print out values of the current jobs SLURM environment variables
  
scontrol show job ${SLURM_JOB_ID}                              # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch

# Note: Run dos2unix {filename} if sbatch DOS line break error occurs
