# Repositories
TOP = ../../Submodules/
AMREX_HOME         = ${TOP}/amrex
IAMR_HOME          = ${TOP}/IAMR
PELELM_HOME        = ${TOP}/PeleLM
PELE_PHYSICS_HOME  = ${TOP}/PelePhysics
AMREX_HYDRO_HOME   = ${TOP}/AMReX-Hydro

# AMREX options
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compiler / parrallel paradigms
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# CVODE
USE_SUNDIALS_PP = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = drm19

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_Model       = Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_Model = Simple

Blocs   := .

CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
