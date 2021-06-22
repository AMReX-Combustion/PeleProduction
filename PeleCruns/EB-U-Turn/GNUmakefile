TOP := $(shell pwd)/../..
PELEC_HOME ?= $(TOP)/Submodules/PeleC
PELE_PHYSICS_HOME ?= $(TOP)/Submodules/PelePhysics
AMREX_HOME ?= $(TOP)/Submodules/amrex
AMREX_HYDRO_HOME ?= $(TOP)/Submodules/AMReX-Hydro

# AMReX
DIM = 3
COMP = clang
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_DPCPP = FALSE

# Debugging
DEBUG = TRUE
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
USE_REACT = FALSE
USE_FORCING = FALSE
USE_EB = TRUE
Eos_dir := Fuego
#Transport_dir := EGLib
Transport_dir := Simple
Chemistry_Model := LiDryer

# GNU Make
Bpack := ./Make.package TurbInflow/Make.package
Blocs := . TurbInflow
include $(PELEC_HOME)/Exec/Make.PeleC
