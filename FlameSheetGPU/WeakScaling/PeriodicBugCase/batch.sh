#!/bin/bash
#BSUB -P cmb138
#BSUB -W 0:30
#BSUB -nnodes 64
#BSUB -alloc_flags gpumps
#BSUB -J FlameSheet-PerBug-0064

module load cmake gcc cuda
jsrun -r 6 -a 1 -g 1 -c 7 -l GPU-CPU -d packed -b rs --smpiargs="-gpu" ./PeleLM3d.gnu.TPROF.MPI.CUDA.ex inputs.3d max_step=2 2>&1 | tee log

# n: # resources sets, where in each resource set:
# a: # ranks
# c: # cores
# g: # gpus
