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
USE_WBAR        = FALSE
VERBOSE         = TRUE
PROFILE         = FALSE
COMP            = gnu
#COMP            = intel
#FCOMP           = intel
COMM_PROFILE    = FALSE

# Enable flctin namelist
USE_FLCTS       = TRUE

#Chemistry_Model = DODECANE_LU_RED
Chemistry_Model = dodecane_wang
#fEXE_sources += getrates.f
fEXE_sources += wang_dodecane_phi_1000.f

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

