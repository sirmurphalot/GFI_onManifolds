#!/bin/bash

#SBATCH --job-name=ar1
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30g
#SBATCH -t 4-00:00:00

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r CMH_ar1 -logfile mycode.out
