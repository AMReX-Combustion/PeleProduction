# Repositories
TOP = $(shell pwd)/../../Submodules
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
VERBOSE         = TRUE
TINY_PROFILE    = TRUE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD = FALSE
USE_EB = FALSE

# PelePhysics
Chemistry_Model = multi_dechep
Eos_Model	= Fuego
Transport_Model = Simple

USE_PARTICLES = TRUE
SPRAY_FUEL_NUM = 2

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
