#!/bin/bash

#SBATCH --time=00:30:00 
#SBATCH --array=1-100

module load R
Rscript standard-bootstrap-comparison.R 500 $SLURM_ARRAY_TASK_ID
