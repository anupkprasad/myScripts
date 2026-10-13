#!/bin/bash
#$ -N L40k998
#$ -q gpuq
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -o test$JOB_ID.out
#$ -e test$JOB_ID.err

mpirun -np 1 /share/apps/namd_gpu/Linux-x86_64-g++/namd2 +p12 +idlepoll prod.inp > prod.log
