TOP = /gpfs/alpine/scratch/lesclapez/cmb138/PeleProduction/Submodules
AMREX_HOME = ${TOP}/amrex
IAMR_HOME = ${TOP}/IAMR
PELELM_HOME = ${TOP}/PeleLM
PELE_PHYSICS_HOME = ${TOP}/PelePhysics
#
# Build configuration
#

DIM             = 3
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE

USE_CUDA        = FALSE

TINY_PROFILE = TRUE

USE_EB = TRUE

USE_SUNDIALS_PP = TRUE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  # use KLU sparse features -- only useful if CVODE is used on CPU
  USE_KLU_PP = FALSE
  USE_CUDA_SUNDIALS_PP=FALSE
endif

# Choose model (from list below), and pmf file
Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Constant

Blocs   := .

CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
