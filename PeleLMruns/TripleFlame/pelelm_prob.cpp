
#include <pelelm_prob.H>
#include <mechanism.h>

namespace ProbParm
{
    AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real splitx = - 0.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real midtanh = 0.001;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real widthtanh = 0.001;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real Zst = 0.055;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real T_in = 300.0;
    AMREX_GPU_DEVICE_MANAGED  amrex::Real V_in = 0.4;

    AMREX_GPU_DEVICE_MANAGED int fuelID = -1;
    AMREX_GPU_DEVICE_MANAGED int bathID = -1;
    AMREX_GPU_DEVICE_MANAGED int oxidID = -1;
} // namespace ProbParm

extern "C" {
    void amrex_probinit (const int* init,
                         const int* name,
                         const int* namelen,
                         const amrex_real* problo,
                         const amrex_real* probhi)
    {
        amrex::ParmParse pp("prob");

        pp.query("P_mean", ProbParm::P_mean);
        pp.query("Zst",    ProbParm::Zst);
        pp.query("T_in",   ProbParm::T_in);
        pp.query("V_in",   ProbParm::V_in);

        ProbParm::splitx = 0.5 * (problo[0] + probhi[0]);
        ProbParm::midtanh = 0.6 * (problo[0] + probhi[0]);
        ProbParm::widthtanh = 0.05 * (problo[0] + probhi[0]); 

        // TODO: somewhat hard coded bath, fuel and oxid IDs 
        // should exist somewhere in PeleLM.
        ProbParm::bathID = N2_ID;  
        ProbParm::fuelID = CH4_ID;  
        ProbParm::oxidID = O2_ID;  
    }
}
