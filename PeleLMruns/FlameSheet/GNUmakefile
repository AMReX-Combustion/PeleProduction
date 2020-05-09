TOP = ../../Submodules
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics
PELELM_HOME ?= ${TOP}/PeleLM
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics

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
TINY_PROFILE    = FALSE
USE_CUDA        = FALSE
USE_CUDA        = TRUE

ifeq ($(USE_CUDA),TRUE)
  USE_SUNDIALS_PP = TRUE
  USE_CUDA_CVODE = TRUE
  USE_KLU_PP = FALSE
else
  USE_SUNDIALS_PP = FALSE
  USE_CUDA_CVODE = FALSE
  USE_KLU_PP = FALSE
endif

#USE_WBAR = TRUE
#DEFINES+=-DALWAYS_NEW_J

DO_2S_CONVERGENCE=FALSE
ifeq ($(DO_2S_CONVERGENCE),TRUE)
  Chemistry_Model = chem-CH4-2step
else
 Chemistry_Model = drm19
 #Chemistry_Model = LiDryer
 #Chemistry_Model = LUDME
 #Chemistry_Model = GRI30
endif

# This sets the EOS directory in $(PELE_PHYSICS_HOME)/Eos
Eos_dir     := Fuego

# This sets the network directory in $(PELE_PHYSICS_HOME)/Reactions
Reactions_dir := Fuego

# This sets the transport directory in $(PELE_PHYSICS_HOME)/Transport
Transport_dir := Simple

Blocs   := .

# Local sources
CEXE_sources +=
F90EXE_sources += Prob_nd.F90 probdata.F90
CEXE_headers += 
FEXE_headers += 

ifeq ($(DO_2S_CONVERGENCE),TRUE)
  F90EXE_sources += 2step_Phi0.8_refined_CubicInterp.F90
  DEFINES+=-DINTERP_PMF_AS_POINT
else
  ifeq ($(Chemistry_Model),drm19)
    fEXE_sources += drm19Soln_seed_0.50.f
  else
    ifeq ($(Chemistry_Model),LiDryer)
      fEXE_sources += LuDME_0700.f
    else
      ifeq ($(Chemistry_Model),LUDME)
        fEXE_sources += LuDME_0700.f
      else
        ifeq ($(Chemistry_Model),GRI30)
          fEXE_sources += gri30_070.f
        endif
      endif
    endif
  endif
endif

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
