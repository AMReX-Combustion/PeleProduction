TOP = $(shell pwd)/../..
SUBMODS = ${TOP}/Submodules
AMREX_HOME ?= ${SUBMODS}/amrex

DEBUG	= FALSE
DEBUG	= TRUE
DIM	= 3
COMP    = gcc
PRECISION = DOUBLE
USE_MPI   = TRUE
USE_OMP   = FALSE

EBASE     = main

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Base
include $(AMREX_HOME)/Src/Base/Make.package
vpathdir += $(AMREX_HOME)/Src/Base

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
include $(AMREX_HOME)/Tools/C_util/Make.package
vpathdir += $(AMREX_HOME)/Tools/C_util

INCLUDE_LOCATIONS += $(TOP)/PeleCruns/EB-U-Turn/TurbInflow
include $(TOP)/PeleCruns/EB-U-Turn/TurbInflow/Make.package
vpathdir += $(TOP)/PeleCruns/EB-U-Turn/TurbInflow

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
