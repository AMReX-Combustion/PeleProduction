#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>

#include "PeleC.H"
#include "IndexDefines.H"
#include "PelePhysics.H"
#include "Tagging.H"
#include "ProblemDerive.H"
#include "prob_parm.H"
#include "Constants.H"
#include "Geometry.H"

class EBUTurnFlipped
  : public pele::pelec::Geometry::Register<EBUTurnFlipped>
{
public:
  static std::string identifier() { return "u_turn-flipped"; }

  void build(const amrex::Geometry& geom, const int max_coarsening_level);
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  amrex::Real xlo,
  amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector,
  const ProbParmDevice& prob_parm)
{
  if (prob_parm.pmf_do_average) {
    int lo_loside = 0;
    int lo_hiside = 0;
    int hi_loside = 0;
    int hi_hiside = 0;
    if (xlo < prob_parm.d_pmf_X[0]) {
      lo_loside = 0;
      lo_hiside = 0;
    }
    if (xhi < prob_parm.d_pmf_X[0]) {
      hi_loside = 0;
      hi_hiside = 0;
    }
    if (xlo > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      lo_loside = prob_parm.pmf_N - 1;
      lo_hiside = prob_parm.pmf_N - 1;
    }
    if (xhi > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      hi_loside = prob_parm.pmf_N - 1;
      hi_hiside = prob_parm.pmf_N - 1;
    }
    if (lo_loside == 0) {
      for (unsigned int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if ((xlo > prob_parm.d_pmf_X[i]) && (xlo < prob_parm.d_pmf_X[i + 1])) {
          lo_loside = i;
          lo_hiside = i + 1;
        }
      }
    }
    if (hi_loside == 0) {
      for (unsigned int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if ((xhi > prob_parm.d_pmf_X[i]) && (xhi < prob_parm.d_pmf_X[i + 1])) {
          hi_loside = i;
          hi_hiside = i + 1;
        }
      }
    }
    for (unsigned int j = 0; j < prob_parm.pmf_M; j++) {
      amrex::Real x1 = prob_parm.d_pmf_X[lo_loside];
      amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + lo_loside];
      amrex::Real x2 = prob_parm.d_pmf_X[lo_hiside];
      amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + lo_hiside];
      amrex::Real dydx = 0.0;
      if (lo_loside == lo_hiside)
        dydx = 0.0;
      else
        dydx = (y2 - y1) / (x2 - x1);
      amrex::Real ylo = y1 + dydx * (xlo - x1);
      amrex::Real yhi = 0.0;
      if (lo_loside == hi_loside) {
        yhi = y1 + dydx * (xhi - x1);
        y_vector[j] = 0.5 * (ylo + yhi);
      } else {
        amrex::Real sum = (x2 - xlo) * 0.5 * (ylo + y2);
        x1 = prob_parm.d_pmf_X[hi_loside];
        y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hi_loside];
        x2 = prob_parm.d_pmf_X[hi_hiside];
        y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hi_hiside];
        if (hi_loside == hi_hiside)
          dydx = 0.0;
        else
          dydx = (y2 - y1) / (x2 - x1);
        yhi = y1 + dydx * (xhi - x1);
        sum += (xhi - x1) * 0.5 * (yhi + y1);
        for (int k = lo_hiside; k < hi_loside - 1; k++) {
          sum += (prob_parm.d_pmf_X[k + 1] - prob_parm.d_pmf_X[k]) * 0.5 *
                 (prob_parm.d_pmf_Y[prob_parm.pmf_N * j + k] +
                  prob_parm.d_pmf_Y[prob_parm.pmf_N * j + k + 1]);
        }
        y_vector[j] = sum / (xhi - xlo);
      }
    }
  } else {

    amrex::Real xmid = 0.5 * (xlo + xhi);
    int loside = -1;
    int hiside = -1;
 
    if (xmid < prob_parm.d_pmf_X[0]) {
      loside = 0;
      hiside = 0;
    }
    if (xmid > prob_parm.d_pmf_X[prob_parm.pmf_N - 1]) {
      loside = prob_parm.pmf_N - 1;
      hiside = prob_parm.pmf_N - 1;
    }
    if (loside == -1) {
      for (unsigned int i = 0; i < prob_parm.pmf_N - 1; i++) {
        if (
          (xmid >= prob_parm.d_pmf_X[i]) &&
          (xmid <= prob_parm.d_pmf_X[i + 1])) {
          loside = i;
          hiside = i + 1;
        }
      }
    }
    for (unsigned int j = 0; j < prob_parm.pmf_M; j++) {
      const amrex::Real x1 = prob_parm.d_pmf_X[loside];
      const amrex::Real y1 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + loside];
      const amrex::Real x2 = prob_parm.d_pmf_X[hiside];
      const amrex::Real y2 = prob_parm.d_pmf_Y[prob_parm.pmf_N * j + hiside];
      amrex::Real dydx = 0.0;
      if (loside == hiside)
        dydx = 0.0;
      else
        dydx = (y2 - y1) / (x2 - x1);
      y_vector[j] = y1 + dydx * (xlo - x1);
    }
  }

  
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_initdata(
	    int i,
	    int j,
	    int k,
	    amrex::Array4<amrex::Real> const& state,
	    amrex::GeometryData const& geomdata,
	    ProbParmDevice const& prob_parm)
{  

  if (prob_parm.phi_in < 0.0)
    {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
      const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
      amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
      amrex::Real pert = 0.0;


      if (prob_parm.pertmag > 0.0)
	{
	  pert =
	    prob_parm.pertmag *
	    (1.0 * std::sin(2 * constants::PI() * 4 * x / prob_parm.L[0]) *
	     std::sin(2 * constants::PI() * 5 * y / prob_parm.L[1]) +
	     1.023 *
	     std::sin(2 * constants::PI() * 2 * (x - 0.4598) / prob_parm.L[0]) *
	     std::sin(2 * constants::PI() * 4 * (y - 0.53765) / prob_parm.L[1]) +
	     0.945 *
	     std::sin(2 * constants::PI() * 3 * (x - 0.712435) / prob_parm.L[0]) *
	     std::sin(2 * constants::PI() * 3 * (y - 2.137) / prob_parm.L[1]) +
	     1.017 *
	     std::sin(2 * constants::PI() * 5 * (x - 0.33) / prob_parm.L[0]) *
	     std::sin(2 * constants::PI() * 6 * (y - 1.8) / prob_parm.L[1]) +
	     0.982 *
	     std::sin(2 * constants::PI() * 5 * (x - 1.4234) / prob_parm.L[0]));
	}
      pmf(y + pert, y + pert, pmf_vals, prob_parm);
      const amrex::Real T = pmf_vals[0];
      //const amrex::Real T = prob_parm.T_in;
      const amrex::Real pres = prob_parm.pamb;
      amrex::Real u[3] = {0.0};
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      amrex::Real rho = 0.0;
      amrex::Real energy = 0.0;

      for (int n = 0; n < NUM_SPECIES; n++)
	massfrac[n] = prob_parm.fuel_state[UFS + n - 1]/prob_parm.fuel_state[URHO];

      u[0] = - pmf_vals[1]*tanh(-10.0*(y-0.5*prob_hi[1]+0.1));
      
      auto eos = pele::physics::PhysicsType::eos();
      //eos.X2Y(molefrac, massfrac);
      eos.PYT2RE(pres, massfrac, T, rho, energy);

      state(i, j, k, URHO) = rho;
      AMREX_D_TERM(state(i, j, k, UMX) = rho * u[0];,
                   state(i, j, k, UMY) = rho * u[1];,
                   state(i, j, k, UMZ) = rho * u[2];)
      state(i, j, k, UEINT) = rho * energy;
      state(i, j, k, UEDEN) =
	rho * (energy + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
      state(i, j, k, UTEMP) = T;
      for (int n = 0; n < NUM_SPECIES; n++)
	state(i, j, k, UFS + n) = rho * massfrac[n];

    }
  else
    {
      for (int n = 0; n < NVAR; n++)
	state(i, j, k, n) = prob_parm.fuel_state[n];
    }
}

void pc_prob_close();

using ProblemTags = EmptyProbTagStruct;
using ProblemDerives = EmptyProbDeriveStruct;

#endif
