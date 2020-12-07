#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm
{
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real splitx;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real midtanh;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real widthtanh;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real Zst;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real T_in;
extern AMREX_GPU_DEVICE_MANAGED  amrex::Real V_in;

extern AMREX_GPU_DEVICE_MANAGED int bathID;
extern AMREX_GPU_DEVICE_MANAGED int fuelID;
extern AMREX_GPU_DEVICE_MANAGED int oxidID;
}

#endif
