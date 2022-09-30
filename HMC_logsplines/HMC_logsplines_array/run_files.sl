#!/bin/bash

#SBATCH --job-name=HMC_array
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=1-100
#SBATCH --mem=10g
#SBATCH -t 6-00:00:00
#SBATCH --output=./SLURMOUT/slurm_log_%A-%a.out

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r HMC_logsplines_array -logfile ./logfiles/testingArray_$SLURM_ARRAY_TASK_ID.out
