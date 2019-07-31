PELELM_HOME ?= ../../Submodules/PeleLM
PELE_PHYSICS_HOME ?= ../../Submodules/PelePhysics
IAMR_HOME ?= ../../Submodules/IAMR
AMREX_HOME ?= ../../Submodules/amrex

#
# User-settable variables.
#
PRECISION         = DOUBLE
DEBUG	          = FALSE
DIM    	          = 3
COMP              = intel
USE_OMP           = FALSE
USE_MPI           = TRUE
VERBOSE           = FALSE
NESAP_BENCHMARK   = TRUE

Chemistry_Model = CHEMH

ifeq ($(NESAP_BENCHMARK), TRUE)
  Chemistry_Model = FFCM1_Red
  COMP = intel
  DIM = 3
  PRECISION = DOUBLE
  DEBUG = FALSE
  USE_OMP = FALSE
  USE_MPI = TRUE
  TINY_PROFILE = TRUE
endif

ifeq ($(Chemistry_Model), CHEMH)
  fEXE_sources += H-flame_0500.f
endif
ifeq ($(Chemistry_Model), FFCM1)
  fEXE_sources += CH4-1_O2-2_N2-2p5_H2-5p833.f
endif

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_F.F90 Prob_$(DIM)d.F90
CEXE_headers += 
FEXE_headers += probdata.H MyProb_F.H

include $(PELELM_HOME)/Tools/Make/Make.PeleLM


