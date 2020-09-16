#!/bin/bash
#BSUB -P cmb138
#BSUB -W 0:30
#BSUB -nnodes 2
#BSUB -alloc_flags gpumps
#BSUB -J FlameSheet-SL-0002

module load cmake gcc cuda
set -x
export OMP_NUM_THREADS=7
# -p X, where X = 6 * number of nodes
jsrun -p 12 -c 7 -r 6 -l cpu-cpu -d packed -b rs ./PeleLM3d.gnu.TPROF.MPI.OMP.ex inputs.3d.sl max_step=2 2>&1 | tee log_SL

# n: # resources sets, where in each resource set:
# a: # ranks
# c: # cores
# g: # gpus
