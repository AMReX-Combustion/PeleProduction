#!/bin/bash
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -q debug
#SBATCH -J DMEtest
#SBATCH -t 00:30:00
#SBATCH -A m3406

#OpenMP settings:
export OMP_NUM_THREADS=8
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
#srun -n 4 -c 8 --cpu_bind=cores ./PeleLM3d.gnu.haswell.MPI.OMP.ex inputs.marc >> runlog
srun -n 32 -c 1 --cpu_bind=cores ./PeleLM3d.gnu.haswell.MPI.ex inputs.marc >> runlogMPI

