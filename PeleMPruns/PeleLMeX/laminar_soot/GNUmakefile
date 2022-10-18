# Repositories
TOP = $(shell pwd)/../../../Submodules
AMREX_HOME         := ${TOP}/amrex
IAMR_HOME          := ${TOP}/IAMR
PELELMEX_HOME      := ${TOP}/PeleLMeX
PELE_PHYSICS_HOME  := ${TOP}/PelePhysics
AMREX_HYDRO_HOME   := ${TOP}/AMReX-Hydro
PELEMP_HOME        := ${TOP}/PeleMP

# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
Chemistry_Model = SootReaction
Eos_Model	= Fuego
Transport_Model = Simple

USE_SOOT = TRUE
NUM_SOOT_MOMENTS = 3
#HIPCC_FLAGS += -mllvm -greedy-regclass-priority-trumps-globalness
include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
