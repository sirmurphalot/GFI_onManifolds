#!/bin/bash

#SBATCH --job-name=HMC
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30g
#SBATCH -t 6-00:00:00

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r CHMC_logspline -logfile mycode.out
