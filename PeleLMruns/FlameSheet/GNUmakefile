# use the full path
TOP               = $(shell pwd)/../..
PELELM_HOME       = $(TOP)/Submodules/PeleLM
PELE_PHYSICS_HOME = $(TOP)/Submodules/PelePhysics
IAMR_HOME         = $(TOP)/Submodules/IAMR
AMREX_HOME        = $(TOP)/Submodules/amrex
AMREX_HYDRO_HOME  = $(TOP)/Submodules/AMReX-Hydro

#
# Build configuration
#

# AMREX options
DIM             = 2

# Compiler / parrallel paradigms
COMP            = gnu
USE_MPI         = FALSE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# MISC options
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = TRUE
TINY_PROFILE    = FALSE

# CVODE
USE_SUNDIALS_PP = TRUE
USE_KLU_PP      = FALSE

# PeleLM options
DO_2S_CONVERGENCE=FALSE

ifeq ($(DO_2S_CONVERGENCE),TRUE)
	Chemistry_Model = chem-CH4-2step
	F90EXE_sources += 2step_Phi0.8_refined_CubicInterp.F90
	DEFINES+=-DINTERP_PMF_AS_POINT
else
	# Enable the restart-from-plotfile option (ht.pltfile)
	#  ... default this off because when on, code expects it to be used
	# USE_NEWMECH    = FALSE

	# Choose model (from list below), and pmf file
	Chemistry_Model = drm19
#	fEXE_sources += drm19Soln_seed_0.50.f

	#Chemistry_Model = LIDRYER
	#fEXE_sources += LiDryer_040_01ATM.f

	#Chemistry_Model = LUDME
	#fEXE_sources += LuDME_0700.f

	#Chemistry_Model = GRI30
	#fEXE_sources += gri30_070.f
endif

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

Blocs   := .

CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

ifeq ($(USE_CUDA),TRUE)
CXXFLAGS += -Xptxas --disable-optimizer-constants
endif
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
