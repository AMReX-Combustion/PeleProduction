PELELM_HOME = ../../../Submodules/PeleLM
PELE_PHYSICS_HOME = ../../../Submodules/PelePhysics
IAMR_HOME = ../../../Submodules/IAMR
AMREX_HOME = ../../../Submodules/amrex

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

#######################
USE_SUNDIALS_PP = FALSE

USE_KLU_PP = FALSE
#######################

Reactions_dir   := Fuego
Transport_dir   := Simple
Eos_dir         := Fuego
Chemistry_Model := grimech30

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_nd.F90 probdata.F90
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
