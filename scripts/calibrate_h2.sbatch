#!/bin/bash
#SBATCH --job-name=pknn_calib
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --array=1-24%12
#SBATCH --cpus-per-task=2
#SBATCH --mem=32gb
#SBATCH --time=04:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/calib_%A_%a.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/calib_%A_%a.err
module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/c$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR
Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/calibrate_h2.R
rm -rf $TMPDIR
