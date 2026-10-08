#!/bin/bash
#SBATCH --job-name=pknn_analyze_sparse
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --array=1-180%25
#SBATCH --cpus-per-task=4
#SBATCH --mem=48gb
#SBATCH --time=20:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/analyze_sparse_%A_%a.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/analyze_sparse_%A_%a.err
module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/analyze_sparse$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR
Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/analyze_sparse.R
rm -rf $TMPDIR
