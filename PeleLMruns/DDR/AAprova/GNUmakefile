PELELM_HOME ?= ../../../Submodules/PeleLM
PELE_PHYSICS_HOME ?= ../../../Submodules/PelePhysics
IAMR_HOME ?= ../../../Submodules/IAMR
AMREX_HOME ?= ../../../Submodules/amrex

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

Reactions_dir   := Fuego
Transport_dir   := Simple
Eos_dir         := Fuego
Chemistry_Model := grimech30

Blocs   := .

F90EXE_sources += Prob_$(DIM)d.F90 user_defined_fcts_$(DIM)d.F90 probdata.F90

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
