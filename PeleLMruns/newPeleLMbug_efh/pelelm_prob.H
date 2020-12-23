#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <PeleLM_parm.H>


#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <EOS.H>
#include <AMReX_CONSTANTS.H>


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 PmfData const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx      = geomdata.CellSize();

  const amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> xx = {{AMREX_D_DECL(prob_lo[0] + (i+0.5)*dx[0],
                                                                        prob_lo[1] + (j+0.5)*dx[1],
                                                                        prob_lo[2] + (k+0.5)*dx[2])}};
  
  const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
  const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
  const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  amrex::Real r = std::sqrt( AMREX_D_TERM(x*x, + y*y, +z*z ));
  AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM>1);
  amrex::Real y1 = ( (prob_parm.FK_rad-r) - prob_parm.standoff - 0.5 * dx[1]) * 100.;
  amrex::Real y2 = ( (prob_parm.FK_rad-r) - prob_parm.standoff + 0.5 * dx[1]) * 100.;
  
  if (prob_parm.pertmag > 0.0) 
  {
#if ( AMREX_SPACEDIM == 2 )
    amrex::Real theta = std::atan(y/x);
    amrex::Real pert = prob_parm.pertmag*(1.000 * sin(4*theta)
                                          + 1.023 * sin(2*(0.25 * 2 * Pi))
                                          + 0.945 * sin(3*(0.1 * 2 * Pi))
                                          + 1.017 * sin(5*(0.7 * 2 * Pi))
                                          + .982 * sin(5*(0.45 * 2 * Pi)));
#else
    pert = pertmag;
#endif
    y1 += pert * 100.;
    y2 += pert * 100.;
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};

  PMF::pmf(pmf_data,y1, y2, pmf_vals);

  state(i,j,k,DEF_Temp) = pmf_vals[0];

  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  for (int n=0; n<NUM_SPECIES; ++n) {
    molefrac[n] = pmf_vals[3+n];
  }
  EOS::X2Y(molefrac, massfrac);

  state(i,j,k,Xvel) = 0;
#if ( AMREX_SPACEDIM == 2 ) 
  state(i,j,k,Yvel) = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 ) 
  state(i,j,k,Yvel) = 0;
  state(i,j,k,Zvel) = pmf_vals[1]*1e-2;
#endif

  amrex::Real P_cgs = prob_parm.P_mean * 10.0;

  amrex::Real rho_cgs;
  EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
  state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

  EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
  state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm, 
  PmfData const *pmf_data)
{
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& dBox,
           const int  dir,
           const int  state_comp,
           const int  nComp)
{
}
#endif
