#!/bin/bash
#SBATCH --job-name=GT        # Job name
#SBATCH --partition=gpu_p              # Partition name
#SBATCH --nodes=1                       # Number of compute nodes
#SBATCH --ntasks-per-node=64                       # Number of MPI tasks
#SBATCH --cpus-per-task=1              # Each MPI task gets 1 CPU 
#SBATCH --mem=20G                        # Memory per node
#SBATCH --time=0-10:00:00               # Time limit
#SBATCH --mail-user=akp66103@uga.edu    # Email
###SBATCH --mail-type=END,FAIL            # Mail events
##SBATCH --gres=gpu:1                   # Uncomment if using GPU-enabled NAMD

# Load NAMD module
module load NAMD/2.14-foss-2023a-mpi

# Run NAMD
mpirun namd2 equilibration.inp > min.log

####charmrun +p64 namd2 min.conf > min.log
##namd2 +p64 min.conf > min.log

### If GPU-enabled     namd2 +p64 +idlepoll min.conf > min.log
