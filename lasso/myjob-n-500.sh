#!/bin/bash

#SBATCH --time=02:00:00 
#SBATCH --array=1-100

module load R
Rscript inference-n-500.R $SLURM_ARRAY_TASK_ID
