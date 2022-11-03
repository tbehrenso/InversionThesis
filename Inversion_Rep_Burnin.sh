#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=8:00:00
#SBATCH --job-name=slimjob
#SBATCH --output=/data/users/tbehrens/InversionThesis/slim_output/slim_%j.o
#SBATCH --error=/data/users/tbehrens/InversionThesis/errors_slim/error_slim_%j.e
#SBATCH --mail-user=thomas.behrens@students.unibe.ch
#SBATCH --mail-type=fail

module load R/latest
module load Development/gcc/9.2.1

s=0.1
m=0.001
mu=1e-5
nrep=100
r=1e-6

sim_type=neutral_2pop
dir_name=${sim_type}_s${s}_m${m}_mu${mu}_r${r}

mkdir Outputs/${dir_name}
mkdir Outputs/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}

for r in $(seq 1 $nrep)
do
	./slim -d rep=$r -d mu=$mu -d s=$s -d m=$m -d r=$r -d "dir_name='$dir_name'" slim_scripts/${sim_type}.slim
done


#Rscript Rtest.R