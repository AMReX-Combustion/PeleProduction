PELELM_HOME ?= ../../Submodules/PeleLM
PELE_PHYSICS_HOME ?= ../../Submodules/PelePhysics
IAMR_HOME ?= ../../Submodules/IAMR
AMREX_HOME ?= ../../Submodules/amrex

#
# Build configuration
#

PRECISION       = DOUBLE
DEBUG           = FALSE
DIM             = 3
USE_OMP         = FALSE
USE_MPI         = TRUE
USE_WBAR        = TRUE
VERBOSE         = TRUE
PROFILE         = FALSE
COMP            = gnu
COMM_PROFILE    = FALSE

# Enable flctin namelist
USE_FLCTS       = TRUE

# 35-species n-dodecane mechanism (LTC/HTC)
# Chemistry_Model = dodecane_lu_red
# 53-species n-dodecane mechanism (LTC/HTC)
Chemistry_Model = dodecane_lu
# 56-species n-dodecane mechanism (HTC only)
#Chemistry_Model = dodecane_wang

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib

Blocs   := .

CEXE_sources +=
F90EXE_sources += Prob_$(DIM)d.F90 probdata.F90
CEXE_headers +=
FEXE_headers +=

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

