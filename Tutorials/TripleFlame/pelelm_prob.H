#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelm_prob_parm.H>
#include <pmf_data.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <EOS.H>
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Phi(ProbParm const& prob_parm,
               const amrex::Real phi_in,
               amrex::Real Y[])
{   
   amrex::Real Xt[NUM_SPECIES] = {0.0};   
   amrex::Real a = 0.0;   
   if ( prob_parm.fuelID == CH4_ID ) {   
      a = 2.0;
   }      
   Xt[prob_parm.oxidID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
   Xt[prob_parm.fuelID] = phi_in * Xt[prob_parm.oxidID] / a;   
   Xt[prob_parm.bathID] = 1.0 - Xt[prob_parm.oxidID] - Xt[prob_parm.fuelID];   

   EOS::X2Y(Xt,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Ksi(ProbParm const& prob_parm,
               const amrex::Real ksi_in,
               amrex::Real Y[])
{   
   amrex::Real phi = ksi_in / ( 1.0 - amrex::min(ksi_in,0.999999) ) * (1.0 - prob_parm.Zst ) / prob_parm.Zst;   

   set_Y_from_Phi(prob_parm,phi,Y);
}

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
    air_Y[prob_parm.bathID] = 0.767;  
    air_Y[prob_parm.oxidID] = 0.233;  

    // Get mixture fraction and corresponding composition
    amrex::Real tanhval = 0.5 * ( 1.0 + std::tanh((x - prob_parm.midtanh) / prob_parm.widthtanh));  
    set_Y_from_Ksi(prob_parm,tanhval,massfrac);

    amrex::Real rad = std::sqrt(  (x - prob_parm.splitx)      * (x - prob_parm.splitx) 
                                + (amrex::min(y,y_lo) - y_lo) * (amrex::min(y,y_lo) - y_lo));
    
    amrex::Real Gaussian_T = std::exp( - (rad*rad) / (2.0 * ( Gauss_T_width * ( 1.0 + (y - y_lo) / 0.03)) 
                                                          * ( Gauss_T_width * ( 1.0 + (y - y_lo) / 0.03)) ));  
    state(i,j,k,DEF_Temp) = prob_parm.T_in + (Gauss_maxT-prob_parm.T_in) * Gaussian_T;

    amrex::Real Gaussian_spec = std::exp( - (rad*rad) / (2.0 * ( Gauss_Spec_width * ( 1.0 + (y - y_lo) / 0.03)) 
                                                          * ( Gauss_Spec_width * ( 1.0 + (y - y_lo) / 0.03)) ));  
    for (int n = 0; n < NUM_SPECIES; ++n) {  
       massfrac[n] = (1.0 - Gaussian_spec) * massfrac[n] + Gaussian_spec * air_Y[n];
    }     

    state(i,j,k,Xvel) = 0.0;
    state(i,j,k,Yvel) = prob_parm.V_in;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

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
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  if (sgn == 1) {

    if (!ac_parm.ctrl_active || time <= 0.0) {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = prob_parm.V_in;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = 0.0;
       s_ext[Zvel] = prob_parm.V_in;
#endif
    } else {
       s_ext[Xvel] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
       s_ext[Yvel] = ac_parm.ctrl_V_in + (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
#elif (AMREX_SPACEDIM == 3)
       s_ext[Yvel] = 0.0;
       s_ext[Zvel] = ac_parm.ctrl_V_in + (time - ac_parm.ctrl_tBase) * ac_parm.ctrl_dV;
#endif
    }

    s_ext[DEF_Temp] = prob_parm.T_in;
      
    amrex::Real tanhval = 0.5 * ( 1.0 + std::tanh((x[0] - prob_parm.midtanh) / prob_parm.widthtanh));  

    set_Y_from_Ksi(prob_parm,tanhval,massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

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

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx      = geomdata.CellSize();
    const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
    const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
    const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

}
#endif
