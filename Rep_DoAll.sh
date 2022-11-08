#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=24:00:00
#SBATCH --job-name=slim_R
#SBATCH --output=/data/users/tbehrens/InversionThesis/slim_output/slim_%j.o
#SBATCH --error=/data/users/tbehrens/InversionThesis/errors_slim/error_slim_%j.e
#SBATCH --mail-user=thomas.behrens@students.unibe.ch
#SBATCH --mail-type=fail
#SBATCH --array=0-3

module load R/latest
module load Development/gcc/9.2.1

s=0.1
m=0.001
mu=1e-5
nrep=100
rec=1e-6

sim_types=("neutral_2pop" "locallyAdapted_2pop" "adaptiveInversion_2pop" "inversionLAA_2pop")

sim_type=${sim_types[$SLURM_ARRAY_TASK_ID]}
dir_name=${sim_type}_s${s}_m${m}_mu${mu}_r${rec}
tempdir=$SCRATCH

mkdir Outputs/${dir_name}
mkdir Outputs/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}

for r in $(seq 1 $nrep)
do
	./slim -d rep=$r -d mu=$mu -d s=$s -d m=$m -d rec=$rec -d "tempdir='$tempdir'" -d "dir_name='$dir_name'" slim_scripts/${sim_type}.slim
done

mkdir Plots/${dir_name}
mkdir Plots/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}

Outputs/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}

for i in {5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}
do
	Rscript R_scripts/Inversion_Rep_General.R ${dir_name} $i
done

