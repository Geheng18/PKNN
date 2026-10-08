#!/bin/bash
#SBATCH --job-name=rd_exp
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --array=1-6
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb
#SBATCH --time=48:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/realdata/shared/logs/exp%a_%A.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/realdata/shared/logs/exp%a_%A.err
module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/rd${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $TMPDIR
Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/rd_analyze.R
rm -rf $TMPDIR
