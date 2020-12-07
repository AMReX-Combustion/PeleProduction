#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <EOS.H>
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Phi(const amrex::Real phi_in,
               amrex::Real Y[])
{   
   amrex::Real Xt[NUM_SPECIES] = {0.0};   
   amrex::Real a = 0.0;   
   if ( ProbParm::fuelID == CH4_ID ) {   
      a = 2.0;
   }      
   Xt[ProbParm::oxidID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
   Xt[ProbParm::fuelID] = phi_in * Xt[ProbParm::oxidID] / a;   
   Xt[ProbParm::bathID] = 1.0 - Xt[ProbParm::oxidID] - Xt[ProbParm::fuelID];   

   EOS::X2Y(Xt,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Ksi(const amrex::Real ksi_in,
               amrex::Real Y[])
{   
   amrex::Real phi = ksi_in / ( 1.0 - amrex::min(ksi_in,0.999999) ) * (1.0 - ProbParm::Zst ) / ProbParm::Zst;   

   set_Y_from_Phi(phi,Y);
}

AMREX_GPU_DEVICE
inline
void
pelelm_initdata_newmech (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata)
{

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = state(i,j,k,DEF_first_spec+n);
    }

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }
}

AMREX_GPU_DEVICE
inline
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    const amrex::Real L_y = prob_hi[1] - prob_lo[1];

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    // Hot gaussian properties
    amrex::Real y_lo = 0.01;
    amrex::Real Gauss_maxT = 1740.0;
    amrex::Real Gauss_T_width = 0.0012;
    amrex::Real Gauss_Spec_width = 0.0020;

    // Define air
    amrex::Real air_Y[NUM_SPECIES] = {0.0};
    air_Y[ProbParm::bathID] = 0.767;  
    air_Y[ProbParm::oxidID] = 0.233;  

    // Get mixture fraction and corresponding composition
    amrex::Real tanhval = 0.5 * ( 1.0 + std::tanh((x - ProbParm::midtanh) / ProbParm::widthtanh));  
    set_Y_from_Ksi(tanhval,massfrac);

    amrex::Real rad = std::sqrt(  (x - ProbParm::splitx)      * (x - ProbParm::splitx) 
                                + (amrex::min(y,y_lo) - y_lo) * (amrex::min(y,y_lo) - y_lo));
    
    amrex::Real Gaussian_T = std::exp( - (rad*rad) / (2.0 * ( Gauss_T_width * ( 1.0 + (y - y_lo) / 0.03)) 
                                                          * ( Gauss_T_width * ( 1.0 + (y - y_lo) / 0.03)) ));  
    state(i,j,k,DEF_Temp) = ProbParm::T_in + (Gauss_maxT-ProbParm::T_in) * Gaussian_T;

    amrex::Real Gaussian_spec = std::exp( - (rad*rad) / (2.0 * ( Gauss_Spec_width * ( 1.0 + (y - y_lo) / 0.03)) 
                                                          * ( Gauss_Spec_width * ( 1.0 + (y - y_lo) / 0.03)) ));  
    for (int n = 0; n < NUM_SPECIES; ++n) {  
       massfrac[n] = (1.0 - Gaussian_spec) * massfrac[n] + Gaussian_spec * air_Y[n];
    }     

    state(i,j,k,Xvel) = 0.0;
    state(i,j,k,Yvel) = ProbParm::V_in;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  if (sgn == 1) {

    if (!ACParm::ctrl_active || time <= 0.0) {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = ProbParm::V_in;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = 0.0;
       s_ext[Zvel] = ProbParm::V_in;
#endif
    } else {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = ACParm::ctrl_V_in + (time - ACParm::ctrl_tBase) * ACParm::ctrl_dV;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = 0.0;
       s_ext[Zvel] = ACParm::ctrl_V_in + (time - ACParm::ctrl_tBase) * ACParm::ctrl_dV;
#endif
    }

    s_ext[DEF_Temp] = ProbParm::T_in;
      
    amrex::Real tanhval = 0.5 * ( 1.0 + std::tanh((x[0] - ProbParm::midtanh) / ProbParm::widthtanh));  

    set_Y_from_Ksi(tanhval,massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = ProbParm::P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, s_ext[DEF_Temp], rho_cgs);
    s_ext[Density] = rho_cgs * 1.0e3;

    EOS::TY2H(s_ext[DEF_Temp], massfrac, RhoH_temp);
    s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
    }
  }

}

AMREX_GPU_DEVICE
inline
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& dBox,
           const int  dir,
           const int  state_comp,
           const int  nComp)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

}
#endif
