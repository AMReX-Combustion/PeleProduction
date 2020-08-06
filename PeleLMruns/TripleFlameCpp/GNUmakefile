TOP = /home/lake/users/amfelden/CCSE/0_CODES/PeleProduction/
PELELM_HOME = $(TOP)/Submodules/PeleLM
PELE_PHYSICS_HOME = $(TOP)/Submodules/PelePhysics
IAMR_HOME = $(TOP)/Submodules/IAMR
AMREX_HOME = $(TOP)/Submodules/amrex

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

TINY_PROFILE = FALSE
#USE_WBAR = TRUE

DEFINES+=-DBL_USE_NEWMECH

#######################
USE_SUNDIALS_PP = TRUE
ifeq ($(USE_SUNDIALS_PP), TRUE)
  # use KLU sparse features -- only useful if CVODE is used on CPU
  USE_KLU_PP = FALSE
  # Activates GPU version
  USE_CUDA_SUNDIALS_PP = FALSE
endif
#######################

# Choose model (from list below), and pmf file
Chemistry_Model = heptane_fc

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

Blocs   := .

F90EXE_sources += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
