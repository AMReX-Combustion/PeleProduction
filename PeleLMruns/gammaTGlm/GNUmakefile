TOP = ../../Submodules
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics
PELELM_HOME ?= ${TOP}/PeleLM
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics
PELEC_MP_HOME ?= ${TOP}/PeleC-MP

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

USE_PARTICLES = TRUE
#SPRAY_COMPONENTS = 8
#SPRAY_FUEL_NUM = 1
#DEFINES += -DSPRAY_COMPONENTS=$(SPRAY_COMPONENTS)
DEFINES += -DSPRAY_FUEL_NUM=1

TINY_PROFILE = FALSE

#######################

# Enable ht.velocity_plotfile option
#USE_VELOCITY    = FALSE

Chemistry_Model = air

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
#Transport_dir := EGLib
Transport_dir := Simple

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_nd.F90 probdata.F90
CEXE_headers += prob_parm.H #prob.H
FEXE_headers += 

ifeq ($(USE_PARTICLES),TRUE)
	CEXE_sources += SprayParticlesInitInsert.cpp
	#CEXE_headers += prob.H
endif

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
