#!/bin/bash

#SBATCH --job-name=ar1_shorter
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=40g
#SBATCH -t 10-00:00:00

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r CMH_ar1 -logfile mycode.out
