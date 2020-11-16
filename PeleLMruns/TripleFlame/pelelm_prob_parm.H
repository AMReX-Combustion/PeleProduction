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

extern AMREX_GPU_DEVICE_MANAGED unsigned int pmf_N;
extern AMREX_GPU_DEVICE_MANAGED unsigned int pmf_M;
extern AMREX_GPU_DEVICE_MANAGED bool pmf_do_average;

extern amrex::Gpu::ManagedVector<amrex::Real>* pmf_X;
extern amrex::Gpu::ManagedVector<amrex::Real>* pmf_Y;

extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_X;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real* d_pmf_Y;

extern std::string pmf_datafile;
extern amrex::Vector<std::string> pmf_names;

extern AMREX_GPU_DEVICE_MANAGED int bathID;
extern AMREX_GPU_DEVICE_MANAGED int fuelID;
extern AMREX_GPU_DEVICE_MANAGED int oxidID;

}

#endif
