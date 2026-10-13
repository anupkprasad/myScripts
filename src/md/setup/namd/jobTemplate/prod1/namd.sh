## Queue it will run in
#PBS -N L20k_n77
#PBS -q gpuq
#PBS -l select=1:ncpus=18:accelerator=True:vntype=cray_compute
#PBS -l walltime=120:00:00
#PBS -l place=pack
#PBS -j oe

module load craype-broadwell
module load craype-accel-nvidia60
module load namd/2.14/gpu-10.2
EXEC= /scratch/apps/namd/2.14/gpu/10.2/CRAY-XC-intel/namd2

cd $PBS_O_WORKDIR

time aprun -n 1 -N 1 -d 18 /scratch/apps/namd/2.14/gpu/10.2/CRAY-XC-intel/namd2 +p18 +idlepoll prod1.conf > prod1.log

