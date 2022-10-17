#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=8:00:00
#SBATCH --job-name=slimjob
#SBATCH --output=/data/users/tbehrens/inv_2022/slim_output/slim_%j.o
#SBATCH --error=/data/users/tbehrens/inv_2022/errors_slim/error_slim_%j.e
#SBATCH --mail-user=thomas.behrens@students.unibe.ch
#SBATCH --mail-type=fail

module load R/latest
module load Development/gcc/9.2.1


s=0.01
m=0.001
mu=1e-6
nrep=200

sim_type=neutral_2pop
dir_name=${sim_type}_s${s}_m${m}_mu${mu}

mkdir Outputs/${dir_name}
mkdir Outputs/${dir_name}/{5000,6000,7000}

for r in $(seq 1 $nrep)
do
	./slim -d rep=$r -d mu=$mu -d s=$s -d m=$m -d "dir_name='$dir_name'" slim_scripts/${sim_type}.slim
done


#Rscript Rtest.R