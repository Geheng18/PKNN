#!/bin/bash
#SBATCH --job-name=pknn_xtra
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --array=1-200%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=128gb
#SBATCH --time=48:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/xtra_%A_%a.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/xtra_%A_%a.err
module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/x$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR
Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/analyze_extra.R
rm -rf $TMPDIR
