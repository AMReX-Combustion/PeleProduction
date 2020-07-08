PELELM_HOME       = ../../Submodules/PeleLM
PELE_PHYSICS_HOME = ../../Submodules/PelePhysics
IAMR_HOME         = ../../Submodules/IAMR
AMREX_HOME        = ../../Submodules/amrex

#
# Build configuration
#

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = FALSE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

TINY_PROFILE = FALSE
#USE_WBAR = TRUE

DEFINES+=-DALWAYS_NEW_J

#######################
USE_SUNDIALS_PP = FALSE

USE_KLU_PP = FALSE

#######################

# Enable ht.velocity_plotfile option
USE_VELOCITY    = TRUE

# Enable the restart-from-plotfile option (ht.pltfile)
#  ... default this off because when on, code expects it to be used
# USE_NEWMECH    = FALSE

# Choose model (from list below), and pmf file
 Chemistry_Model = drm19
 fEXE_sources += drm19-p01atm-T298K-phi1.f

#Chemistry_Model = grimech30
#fEXE_sources += gri30_070.f

#Chemistry_Model = grimech30-noArN
#fEXE_sources += grimech30-noArN-CH4-0_70.f

#Chemistry_Model = LUDME
#fEXE_sources += LuDME_0700.f

#Chemistry_Model = GRI30
#fEXE_sources += gri30_070.f

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := EGLib

Blocs   := .. .

CEXE_sources +=
F90EXE_sources += Prob_nd.F90 probdata.F90
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
# Available chemistry models:
# CanuH2 chem-CH4-2step chem-H dme glar gri Hai H-CW header inert Konnov
# LiDryer Marinov prf_ethanol Roehl sandiego smooke usc
