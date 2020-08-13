
#include <pelelm_prob_parm.H>
#include <AMReX_ParmParse.H>
#include <pmf.H>

using namespace amrex;

namespace ProbParm
{
  AMREX_GPU_DEVICE_MANAGED  amrex::Real P_mean = 101325.0;
  AMREX_GPU_DEVICE_MANAGED  amrex::Real standoff = -0.03501061;
  AMREX_GPU_DEVICE_MANAGED  amrex::Real pertmag = 0.0;
  AMREX_GPU_DEVICE_MANAGED  amrex::Real FK_rad = 3.8e-05;
}

extern "C"
{
  void amrex_probinit (const int* init,
                       const int* name,
                       const int* namelen,
                       const amrex_real* problo,
                       const amrex_real* probhi)
  {
    ParmParse pp("prob");

    pp.query("P_mean", ProbParm::P_mean);
    pp.query("standoff", ProbParm::standoff);
    pp.query("pertmag", ProbParm::pertmag);
    pp.query("FK_rad", ProbParm::FK_rad);

    std::string pmf_datafile;
    pp.get("pmf_datafile", pmf_datafile);

    PMF::pmf_do_average = true;
    PMF::read_pmf(pmf_datafile);
  }
}
