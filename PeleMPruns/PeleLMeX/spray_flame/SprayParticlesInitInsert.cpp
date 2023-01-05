
#include "SprayParticles.H"
#include "SprayInjection.H"
#include "pelelm_prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real time,
  amrex::Real dt,
  int nstep,
  int lev,
  int finest_level,
  ProbParm const& prob_parm)
{
  amrex::ignore_unused(nstep, finest_level);

  if (lev != 0) {
    return false;
  }
  for (int jet = 0; jet < m_sprayJets.size(); ++jet) {
    SprayJet* js = m_sprayJets[jet].get();
    if (js->jet_active(time)) {
      sprayInjection(time, js, dt, lev);
    }
  }

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(
  const bool init_parts, ProbParm const& prob_parm)
{
  amrex::ignore_unused(prob_parm);
  m_sprayJets.resize(2);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<SprayJet>(jet_name, Geom(0));
  jet_name = "jet2";
  m_sprayJets[1] = std::make_unique<SprayJet>(jet_name, Geom(0));
  // Start without any particles
  m_injectVel = m_sprayJets[0]->jet_vel();
  return;
}
