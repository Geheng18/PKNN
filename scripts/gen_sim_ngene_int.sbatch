#!/bin/bash
#SBATCH --job-name=pknn_ngenI
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --array=1-32%16
#SBATCH --cpus-per-task=1
#SBATCH --mem=24gb
#SBATCH --time=08:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/ngenI_%A_%a.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/ngenI_%A_%a.err
module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/ngenI$SLURM_ARRAY_TASK_ID
mkdir -p $TMPDIR
Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/gen_sim_ngene_int.R
rm -rf $TMPDIR
