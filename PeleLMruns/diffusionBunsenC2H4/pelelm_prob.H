#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <PelePhysics.H>

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

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
    auto eos = pele::physics::PhysicsType::eos();

    const amrex::Real L_x = prob_hi[0] - prob_lo[0];
    amrex::Real splitx = prob_lo[0] + 0.5 * L_x;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};

    massfrac_air[N2_ID] = 0.767;
    massfrac_air[O2_ID] = 0.233;
    massfrac_jet[N2_ID] = 0.0;
    massfrac_jet[O2_ID] = 0.0;
    massfrac_jet[C2H4_ID] = 1.0;

    state(i,j,k,DEF_Temp) = 300.0;

    // Let's setup a field of Z:
    //amrex::Real zfrac = 0.0;
    //amrex::Real zStoich = 0.06;  // close enough
    //amrex::Real jetCoreWidth = std::max(prob_parm.slot_width - y * prob_parm.slot_width / ( 4*prob_parm.slot_width ),0.0);
    //if ( std::abs(x) <= jetCoreWidth ) {
    //    zfrac = 1.0;
    //} else {
    //    zfrac = 1.0 - std::min((std::abs(x)-jetCoreWidth)/0.003,1.0);
    //}
    //if ( y > 4.1*prob_parm.slot_width) {
    //    amrex::Real dist = std::abs(y-4.1*prob_parm.slot_width);
    //    zfrac *= std::max( (0.003-dist)/0.003,0.0) ;
    //}
    //for (int n = 0; n < NUM_SPECIES; n++){
    //    massfrac[n] = zfrac * massfrac_jet[n] + ( 1.0 - zfrac) * massfrac_air[n];
    //}
    //// Put some heat near stoichio
    amrex::Real gaussr = std::abs(y-0.015);
    state(i,j,k,DEF_Temp) = 300 + (prob_parm.Tstart-300) * std::exp(-gaussr*gaussr/(0.0055*0.0055));

    state(i,j,k,Xvel) = 0.0;
    state(i,j,k,Yvel) = 0.01;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac_air, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,DEF_Temp), massfrac_air, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac_air[n] * state(i,j,k,Density);
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

  const amrex::Real L_x = prob_hi[0] - prob_lo[0];
  amrex::Real splitx = prob_lo[0] + 0.5 * L_x;
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
  auto eos = pele::physics::PhysicsType::eos();

  massfrac_air[N2_ID] = 0.767;
  massfrac_air[O2_ID] = 0.233;
  massfrac_jet[N2_ID] = 0.0;
  massfrac_jet[O2_ID] = 0.0;
  massfrac_jet[C2H4_ID] = 1.0;

  for (int n = 0; n < DEF_NUM_STATE; n++){
     s_ext[n] = 0.0;
  }

  s_ext[DEF_Temp] = 300.0;

  // Get zone
  int zone = -1; 
  if ( std::abs(x[0]) < 0.005 ) { 
     zone = 3; // jet
  } else {
     zone = 1; // coflow
  }
  amrex::Real P_cgs = prob_parm.P_mean * 10.0;

  amrex::Real pert = 0.3*prob_parm.Vin;

  // Time variation
  pert *= std::sin(2*Pi*time/0.006) * sin(2*Pi*(time-0.00166)/0.00033333);

  // Space shape
  pert *= (  1.00 * std::sin(2*Pi*0.7777*x[0]/prob_parm.slot_width)
           + 1.07 * std::sin(2*Pi*1.3333*(x[0]-.00127)/prob_parm.slot_width) );

  s_ext[Xvel] = 0.0;
  amrex::Real rad = 0.0;
  if ( zone == 3 ) {
     rad = std::abs(x[0]-splitx);
     s_ext[Yvel] = (prob_parm.Vin+pert) * ( 1.0 - (rad/prob_parm.slot_width) * (rad/prob_parm.slot_width) );
     s_ext[DEF_Temp] = 300.0;
  } else if ( zone == 1 ) {
     rad = std::abs(x[0]-splitx);
     s_ext[Yvel] = prob_parm.Vcoflow * ( 1.0 + (rad)/0.025);
     if ( rad <= prob_parm.slot_width*1.4 ) s_ext[Yvel] = 0.0;
  }

  if ( zone == 1 ) {
     amrex::Real rho_cgs = 0.0;
     amrex::Real RhoH_temp = 0.0;
     eos.PYT2R(P_cgs, massfrac_air, s_ext[DEF_Temp], rho_cgs);
     s_ext[Density] = rho_cgs * 1.0e3;

     eos.TY2H(s_ext[DEF_Temp], massfrac_air, RhoH_temp);
     s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

     for (int n = 0; n < NUM_SPECIES; n++) {
       s_ext[DEF_first_spec+n] = massfrac_air[n] * s_ext[Density];
     }
   } else if ( zone == 3 ) {
     amrex::Real rho_cgs = 0.0;
     amrex::Real RhoH_temp = 0.0;
     eos.PYT2R(P_cgs, massfrac_jet, s_ext[DEF_Temp], rho_cgs);
     s_ext[Density] = rho_cgs * 1.0e3;

     eos.TY2H(s_ext[DEF_Temp], massfrac_jet, RhoH_temp);
     s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion

     for (int n = 0; n < NUM_SPECIES; n++) {
       s_ext[DEF_first_spec+n] = massfrac_jet[n] * s_ext[Density];
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
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  nComp)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();
   const amrex::Real* prob_hi = geomdata.ProbHi();
   const amrex::Real L_x = prob_hi[0] - prob_lo[0];
   const amrex::Real splitx = prob_lo[0] + 0.5 * L_x;
   const amrex::Real* dx      = geomdata.CellSize();
   const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];
   const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];
   const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

   // Need to pass the prob_parm in there.
   //amrex::Real slot_width = 0.005;

   //amrex::Real rad = std::abs(x-splitx);
   //if (  rad >= slot_width &&
   //      rad <= slot_width*1.4 &&
   //      j <= dBox.smallEnd(1)) {
   //   for (int n = 0; n < nComp; n++ ) {
   //      beta(i,j,k,n) = 0.0;
   //   }
   //}
}
#endif
