#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
    amrex::Real P_mean   = 101325.0_rt;
    amrex::Real standoff = 0.0_rt;
    amrex::Real pertmag  = 0.0004_rt;
    amrex::Real Vin = 0.8;
    amrex::Real Vcoflow = 0.1;
    amrex::Real slot_width = 0.005;
    amrex::Real Tstart = 1500;
};
#endif
