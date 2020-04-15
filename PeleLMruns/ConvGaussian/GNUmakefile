TOP = ../../Submodules
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELELM_HOME ?= ${TOP}/PeleLM
PELE_PHYSIC_HOME ?= ${TOP}/PelePhysics

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = TRUE
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

USE_WBAR        = FALSE
USE_VELOCITY    = FALSE
USE_NEWMECH     = FALSE

# Choose model (from list below), and pmf file
Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
#Transport_dir := EGLib
Transport_dir := Constant

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_nd.F90 probdata.F90
CEXE_headers += 
FEXE_headers += 

VPATH_LOCATIONS += .
include $(PELELM_HOME)/Tools/Make/Make.PeleLM
