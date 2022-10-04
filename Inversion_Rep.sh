#!/bin/bash

cd /home/tbehrens/tbehr/Desktop/Thesis

slim=/mnt/c/msys64/mingw64/bin/slim.exe

s=0.1
m=0.01

for r in $(seq 1 100)
do
	$slim -d rep=$r -d mu=1e-6 -d s=$s -d m=$m slim_scripts/inversion_2pop_cmd.slim
done

#Rscript Rtest.R test1 test2