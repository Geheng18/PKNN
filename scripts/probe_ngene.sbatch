#!/bin/bash
#SBATCH --job-name=pknn_probe
#SBATCH --account=lucienq
#SBATCH --qos=lucienq-b
#SBATCH --cpus-per-task=4
#SBATCH --mem=120gb
#SBATCH --time=04:00:00
#SBATCH --output=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/probe_%j.out
#SBATCH --error=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/logs/probe_%j.err
module load R/4.1 2>/dev/null
export PATH=/apps/R/4.1/bin:$PATH
export HPC_SINGULARITY_IMAGE=/apps/R/4.1/container
export TMPDIR=/home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/tmp/probe$SLURM_JOB_ID
mkdir -p $TMPDIR
Rscript /home/heng.ge/blue/heng.ge/KNN_VC_Selection_Revision/scripts/probe_ngene.R
rm -rf $TMPDIR
