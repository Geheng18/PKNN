#!/bin/bash
#SBATCH --job-name=pknn_bench
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --array=1-28%6
#SBATCH --cpus-per-task=4
#SBATCH --mem=160gb
#SBATCH --time=24:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/bench_%A_%a.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/bench_%A_%a.err

module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/b$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR

Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/benchmark_NP.R
rm -rf $TMPDIR
