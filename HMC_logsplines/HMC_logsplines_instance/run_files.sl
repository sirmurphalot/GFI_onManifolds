#!/bin/bash

#SBATCH --job-name=HMC_picture
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4g
#SBATCH -t 02-00:00:00
#SBATCH --output=./SLURMOUT/slurm_log_%A-%a.out

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r CHMC_logspline -logfile ./logfiles/mycode.out
