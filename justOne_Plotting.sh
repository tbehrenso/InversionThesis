#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=13:00:00
#SBATCH --job-name=R_Fst_Weighted
#SBATCH --output=/data/users/tbehrens/InversionThesis/slim_output/RPlot_%j.o
#SBATCH --error=/data/users/tbehrens/InversionThesis/errors_slim/error_RPlot_%j.e
#SBATCH --mail-user=thomas.behrens@students.unibe.ch
#SBATCH --mail-type=end,fail

module load R/latest

s=0.1
m=0.01
mu=1e-5
nrep=100
r=1e-6
rec=1e-6

sim_type=inversionLAA_2pop
dir_name=${sim_type}_s${s}_m${m}_mu${mu}_r${rec}

Outputs/${dir_name}/{10000,15000}

#mkdir Plots/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}
#Outputs/${dir_name}/{5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000}

for i in {5000,10000,15000}
do
	#Rscript R_scripts/justFstHaplo_Hudson.R ${dir_name} $i
  Rscript R_scripts/justFstHaplo.R ${dir_name} $i
done
