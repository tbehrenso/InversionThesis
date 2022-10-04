#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00
#SBATCH --job-name=slimjob
#SBATCH --error=/data/users/tbehrens/error_slim_%j.e
#SBATCH --mail-user=thomas.behrens@students.unibe.ch
#SBATCH --mail-type=fail

module load R/latest

s=0.1
m=0.01

for r in $(seq 1 100)
do
	./slim -d rep=$r -d mu=1e-6 -d s=$s -d m=$m slim_scripts/inversion_2pop_cmd.slim
done


#Rscript --no-save Rtest.R