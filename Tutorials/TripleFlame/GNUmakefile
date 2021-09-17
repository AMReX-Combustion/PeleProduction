TOP = $(shell pwd)/../..
PELELM_HOME       ?= $(TOP)/Submodules/PeleLM
PELE_PHYSICS_HOME ?= $(TOP)/Submodules/PelePhysics
IAMR_HOME         ?= $(TOP)/Submodules/IAMR
AMREX_HOME        ?= $(TOP)/Submodules/amrex
AMREX_HYDRO_HOME  ?= $(TOP)/Submodules/AMReX-Hydro

#
# Build configuration
#
DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

TINY_PROFILE = FALSE

#######################
# Choose model (from list below), and pmf file
Chemistry_Model = drm19

# This sets the EOS model in $(PELE_PHYSICS_HOME)/Eos
Eos_Model      := Fuego

# This sets the transport model in $(PELE_PHYSICS_HOME)/Transport
Transport_Model := Simple

Blocs   := .

F90EXE_sources += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
