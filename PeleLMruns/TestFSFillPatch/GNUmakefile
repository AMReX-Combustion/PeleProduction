# Repositories
TOP = $(PWD)/../../Submodules/
AMREX_HOME         = ${TOP}/amrex
IAMR_HOME          = ${TOP}/IAMR
PELELM_HOME        = ${TOP}/PeleLM
PELE_PHYSICS_HOME  = ${TOP}/PelePhysics
AMREX_HYDRO_HOME   = ${TOP}/AMReX-Hydro

#
# Build configuration
#

# AMREX options
DIM             = 3

# Compiler / parrallel paradigms
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = TRUE
USE_HIP         = FALSE

# MISC options
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE

# This is the chemistry model
Chemistry_Model = dodecane_lu

# This sets the EOS model to one in $(PELE_PHYSICS_HOME)/Eos
Eos_Model       = Fuego

# This sets the transport to one in $(PELE_PHYSICS_HOME)/Transport
Transport_Model = Simple

Blocs   := .
CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
