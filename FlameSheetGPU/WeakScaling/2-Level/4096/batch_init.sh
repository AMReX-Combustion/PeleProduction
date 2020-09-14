#!/bin/bash
#BSUB -P cmb138
#BSUB -W 0:30
#BSUB -nnodes 4096
#BSUB -alloc_flags gpumps
#BSUB -J FlameSheet-init-4096

module load cmake gcc cuda
jsrun -r 6 -a 1 -g 1 -c 7 -l GPU-CPU -d packed -b rs --smpiargs="-gpu" ./PeleLM3d.gnu.TPROF.MPI.CUDA.ex inputs.3d max_step=0 2>&1 | tee log_init

# n: # resources sets, where in each resource set:
# a: # ranks
# c: # cores
# g: # gpus
