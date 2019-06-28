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
COMP              = gnu
USE_OMP           = FALSE
USE_MPI           = TRUE
VERBOSE           = FALSE

Eos_dir           = Fuego
Reactions_dir     = Fuego
Transport_dir     = EGLib

#Chemistry_Model = CHEMH
#fEXE_sources += H-flame_0500.f

Chemistry_Model = FFCM1_Red
fEXE_sources += CH4-1_O2-2_N2-2p5_H2-5p833.f

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_$(DIM)d.F90 probdata.F90
CEXE_headers += 
FEXE_headers +=

include $(PELELM_HOME)/Tools/Make/Make.PeleLM


