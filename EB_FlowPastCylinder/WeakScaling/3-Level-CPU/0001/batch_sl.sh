#!/bin/bash
#BSUB -P cmb138
#BSUB -W 0:30
#BSUB -nnodes 1
#BSUB -alloc_flags gpumps
#BSUB -J FlowPastCyl-3L-0001

module load cmake gcc cuda
set -x
##export OMP_NUM_THREADS=7
##jsrun -p 6 -c 7 -r 6 -l cpu-cpu -d packed -b rs ./PeleLM3d.gnu.TPROF.MPI.OMP.ex inputs.3d max_step=2 2>&1 | tee log_3L
export OMP_NUM_THREADS=1
jsrun -p 6 -c 1 -r 6 -l cpu-cpu -d packed -b rs ./PeleLM3d.gnu.TPROF.MPI.OMP.ex inputs.3d max_step=2 2>&1 | tee log_3L

# n: # resources sets, where in each resource set:
# a: # ranks
# c: # cores
# g: # gpus
