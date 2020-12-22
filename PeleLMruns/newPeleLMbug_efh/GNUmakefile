TOP = ../../Submodules
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = TRUE
USE_CUDA        = FALSE
USE_SUNDIALS_PP = TRUE
TINY_PROFILE    = FALSE

USE_VELOCITY    = TRUE

#USE_WBAR = TRUE
#DEFINES+=-DALWAYS_NEW_J

ifeq ($(USE_CUDA), TRUE)
  USE_KLU_PP = FALSE
endif


Chemistry_Model = drm19
Chemistry_Model = LiDryer
Eos_dir         = Fuego
Reactions_dir   = Fuego
Transport_dir   = Simple

include $(PELE_PHYSICS_HOME)/Utility/Make.package
include $(PELELM_HOME)/Tools/Make/Make.PeleLM
