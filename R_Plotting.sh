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

s=0.01
m=0.001
mu=1e-6
nrep=100

sim_type=inversionLAA_2pop
dir_name=${sim_type}_s${s}_m${m}_mu${mu}

Outputs/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}

for i in {5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}
do
	....
done


#Rscript Rtest.R