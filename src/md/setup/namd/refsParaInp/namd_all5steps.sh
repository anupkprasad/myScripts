#!/bin/bash
#$ -N Dk_440
#$ -q gpuq
#$ -S /bin/bash
#$ -cwd
#$ -pe smp 12
#$ -o test$JOB_ID.out
#$ -e test$JOB_ID.err

equi_prefix=step6.%d_equilibration
prod_prefix=step7_production
prod_step=step7
#Running equilibration steps
cnt=1
cntmax=6

while [ $cnt -le ${cntmax} ]; do
    step=`printf ${equi_prefix} ${cnt}`
    mpirun -np 1 /share/apps/namd_gpu/Linux-x86_64-g++/namd2 +p12 +idlepoll ${step}.inp > ${step}.out
    echo $cnt
    ((cnt++))
done
