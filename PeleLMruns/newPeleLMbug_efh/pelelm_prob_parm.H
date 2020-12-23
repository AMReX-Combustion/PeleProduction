#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real standoff = -0.03501061_rt;
  amrex::Real pertmag = 0.0_rt;
  amrex::Real FK_rad = 3.8e-05_rt;
};
#endif
