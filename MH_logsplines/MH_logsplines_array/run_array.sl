#!/bin/bash

#SBATCH --job-name=logArray
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=3g
#SBATCH --array=1-100
#SBATCH -t 1-00:00:00
#SBATCH --output=./SLURMOUT/slurm_log_%A-%a.out

module add matlab
matlab -nodesktop -nosplash -singleCompThread -r CMH_logsplines_array -logfile ./Rlogfiles/testingArray_$SLURM_ARRAY_TASK_ID.Rout
