#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=1:00:00
#SBATCH --job-name=slimjob
#SBATCH --output=/data/users/tbehrens/InversionThesis/slim_output/slim_%j.o
#SBATCH --error=/data/users/tbehrens/InversionThesis/errors_slim/error_slim_%j.e
#SBATCH --mail-user=thomas.behrens@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --array=0-3

sim_types=("neutral_2pop" "locallyAdapted_2pop" "adaptiveInversion_2pop" "inversionLAA_2pop")

sim_type=${sim_types[$SLURM_ARRAY_TASK_ID]}
dir_name=${sim_type}
tempdir=$SCRATCH


mkdir ${dir_name}

