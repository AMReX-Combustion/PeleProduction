#include <PeleLM.H>
#include <pelelm_prob.H>
#include <PMF.H>

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T_gas", PeleLM::prob_parm->T_init);
}
}
