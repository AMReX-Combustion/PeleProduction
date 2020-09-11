
#include <SprayParticles.H>
#include <AMReX_Particles.H>
//#include <PeleC.H>
#include <PeleLM.H>
//#include "prob.H"

using namespace amrex;

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p0 = 1.013e6; // [erg cm^-3]
AMREX_GPU_DEVICE_MANAGED amrex::Real T0 = 300.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real rho0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real v0 = 0.;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_vel = 2800.;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_dia = .00017; // converted from 1e-2 to 1e-4 for PeleLM
AMREX_GPU_DEVICE_MANAGED amrex::Real part_mean_dia = 5E-4; 
AMREX_GPU_DEVICE_MANAGED amrex::Real part_stdev_dia = 1e-4;
AMREX_GPU_DEVICE_MANAGED amrex::Real mass_flow_rate = 2.0349e-2; //2.3;
AMREX_GPU_DEVICE_MANAGED unsigned int inject_N = 0;
amrex::Gpu::ManagedVector<amrex::Real>* inject_time = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* inject_mass = nullptr;
amrex::Gpu::ManagedVector<amrex::Real>* inject_vel = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_time = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_mass = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real* d_inject_vel = nullptr;
AMREX_GPU_DEVICE_MANAGED amrex::Real part_rho = 0.693;
AMREX_GPU_DEVICE_MANAGED amrex::Real part_temp = 300.;
  //AMREX_GPU_DEVICE_MANAGED amrex::Real Y_O2 = 0.233;
  //AMREX_GPU_DEVICE_MANAGED amrex::Real Y_N2 = 0.767;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_start_time = 0;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_end_time = 10000.;
AMREX_GPU_DEVICE_MANAGED amrex::Real spray_angle = 5.;
AMREX_GPU_DEVICE_MANAGED amrex::Real jet_cent[AMREX_SPACEDIM] = {0.0};
} // namespace ProbParm  

int
interpolateInjectTime(const Real& time)
{
  const int nvals = ProbParm::inject_N;
  int i = 0;
  Real ctime = ProbParm::d_inject_time[i];
  while (ctime < time) {
    ctime = ProbParm::d_inject_time[++i];
  }
  return i;
}

bool
SprayParticleContainer::insertParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  if (lev != 0) return false;
  if (time < ProbParm::jet_start_time
      || time > ProbParm::jet_end_time) return false;

  Print() << "Injecting particles \n";

  // Number of particles per parcel
  const Real num_ppp = 1;//m_parcelSize;
  const Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSize();
  RealVect dom_len(AMREX_D_DECL(geom.ProbLength(0),
                                geom.ProbLength(1),
                                geom.ProbLength(2)));
  Real mass_flow_rate = ProbParm::mass_flow_rate;
  Real jet_vel = ProbParm::jet_vel;
  Real jet_dia = ProbParm::jet_dia;
  Real jr2 = jet_dia*jet_dia/4.; // Jet radius squared
#if AMREX_SPACEDIM == 3
  Real jet_area = M_PI*jr2;
#else
  Real jet_area = jet_dia;
#endif
  Real part_temp = ProbParm::part_temp;
  Real part_rho = ProbParm::part_rho;
  if (ProbParm::inject_N > 0) {
    const int time_indx = interpolateInjectTime(time);
    const Real time1 = ProbParm::d_inject_time[time_indx - 1];
    const Real time2 = ProbParm::d_inject_time[time_indx];
    const Real mf1 = ProbParm::d_inject_mass[time_indx - 1];
    const Real mf2 = ProbParm::d_inject_mass[time_indx];
    const Real invt = (time - time1)/(time2 - time1);
    mass_flow_rate = mf1 + (mf2 - mf1)*invt;
    const Real jv1 = ProbParm::d_inject_vel[time_indx - 1];
    const Real jv2 = ProbParm::d_inject_vel[time_indx];
    jet_vel = jv1 + (jv2 - jv1)*invt;
  }
  //if (jet_vel*dt/dx[0] > 0.5) {
  //  Real max_vel = dx[0]*0.5/dt;
  //  std::string warn_msg = "Injection velocity of " + std::to_string(jet_vel) 
  //    + " is reduced to maximum " + std::to_string(max_vel);
  //  Real m_injectVel = jet_vel;
  //  jet_vel = max_vel;
  //  amrex::Warning(warn_msg);
  //}
  Real part_dia = ProbParm::part_mean_dia;
  Real part_stdev = ProbParm::part_stdev_dia;
  Real stdsq = part_stdev*part_stdev;
  Real meansq = part_dia*part_dia;
  Real log_mean = 2.*std::log(part_dia) - 0.5*std::log(stdsq + meansq);
  Real log_stdev = std::sqrt(amrex::max(-2.*std::log(part_dia)
                                        + std::log(stdsq + meansq), 0.));
  Real Pi_six = M_PI/6.;
  Real spray_angle = ProbParm::spray_angle;
  Real lo_angle = -0.5*spray_angle;

  Print() << "Injecting particles, entering mfiter \n";

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const RealBox& temp = RealBox(bx,geom.CellSize(),geom.ProbLo());
    const Real* xlo = temp.lo();
    const Real* xhi = temp.hi();
    RealVect box_len(AMREX_D_DECL(temp.length(0), temp.length(1), 0.));
    Gpu::HostVector<ParticleType> host_particles;
#ifdef USE_SPRAY_SOA
    std::array<Gpu::HostVector<Real>, NAR_SPR > host_real_attribs;
#endif
    if (xlo[2] == plo[2]) {
      // Box locations relative to jet center
      const RealVect xloJ(AMREX_D_DECL(xlo[0] - ProbParm::jet_cent[0], 
                                       xlo[1] - ProbParm::jet_cent[1],
				       plo[2]));
      const RealVect xhiJ(AMREX_D_DECL(xhi[0] - ProbParm::jet_cent[0],
                                       xhi[1] - ProbParm::jet_cent[1],
				       plo[2]));
      Real cur_jet_area = 0.;
      Real testdx = dx[0]/100.;
      Real testdx2 = testdx*testdx;
      Real curx = xloJ[0];
      // Loop over each cell and check how much overlap there is with the jet
#if AMREX_SPACEDIM == 3
      while (curx < xhiJ[0]) {
        Real curz = xloJ[1];
        while (curz < xhiJ[1]) {
          Real r2 = curx*curx + curz*curz;
          if (r2 <= jr2) cur_jet_area += testdx2;
          curz += testdx;
        }
        curx += testdx;
      }
#else
      while (curx < xhiJ[0]) {
        Real r2 = curx*curx;
        if (r2 <= jr2) cur_jet_area += testdx;
        curx += testdx;
      }
#endif
      Real jet_perc = cur_jet_area/jet_area;
      Real perc_mass = jet_perc*mass_flow_rate*dt;
      Real total_mass = 0.;
      while (total_mass < perc_mass) {
        RealVect part_loc(AMREX_D_DECL(xlo[0] + amrex::Random()*box_len[0],
                                       xlo[1] + amrex::Random()*box_len[1],
				       plo[2]));
	//Print() << "Injecting particles " << part_loc[0] << " " << part_loc[1] << " " << part_loc[2] << " \n";
	//Print() << "Injecting particles " << box_len[0] << " " << box_len[1] << " " << box_len[2] << " \n";
	//Print() << "Injecting particles " << total_mass << " \n";
        Real r2 = AMREX_D_TERM(std::pow(part_loc[0] - ProbParm::jet_cent[0], 2),,
                               + std::pow(part_loc[1] - ProbParm::jet_cent[1], 2));
	//Print() << "Injecting particles " << r2 << " " << jr2 << " \n";
        if (r2 <= jr2) {
	  //Print() << "Actually injecting particle \n";

          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          Real theta = lo_angle + spray_angle*amrex::Random();
#if AMREX_SPACEDIM == 3
          Real theta2 = 2.*M_PI*amrex::Random();
#else
          Real theta2 = 0.;
#endif
          //Real x_vel = jet_vel*std::sin(theta)*std::cos(theta2);
          //Real y_vel = jet_vel*std::cos(theta);
          //Real z_vel = jet_vel*std::sin(theta)*std::sin(theta2);
	  Real x_vel = 2.*(amrex::Random()-.5)*jet_vel*.05;
	  Real y_vel = 2.*(amrex::Random()-.5)*jet_vel*.05;
	  Real z_vel = jet_vel * (1. + 2.*(amrex::Random()-.5)*.05);
	  //Print() << "Injecting particle  z_vel " << z_vel <<" \n";
#ifdef USE_SPRAY_SOA
          AMREX_D_TERM(host_real_attribs[PeleLM::pstateVel].push_back(x_vel);,
                       host_real_attribs[PeleLM::pstateVel+1].push_back(y_vel);,
                       host_real_attribs[PeleLM::pstateVel+2].push_back(z_vel););
#else
          AMREX_D_TERM(p.rdata(PeleLM::pstateVel) = x_vel;,
                       p.rdata(PeleLM::pstateVel+1) = y_vel;,
                       p.rdata(PeleLM::pstateVel+2) = z_vel;);
#endif
          Real cur_dia = amrex::RandomNormal(log_mean, log_stdev);
          // Use a log normal distribution
          cur_dia = std::exp(cur_dia);
          Real part_y = plo[2] + amrex::Random()*dt*y_vel;
          AMREX_D_TERM(p.pos(0) = part_loc[0];,
                       p.pos(1) = part_loc[1];,
                       p.pos(2) = part_y;);

#ifdef USE_SPRAY_SOA
          host_real_attribs[PeleLM::pstateT].push_back(part_temp);
          host_real_attribs[PeleLM::pstateDia].push_back(cur_dia);
          host_real_attribs[PeleLM::pstateRho].push_back(part_rho);
          host_real_attribs[PeleLM::pstateY].push_back(1.);
          for (int sp = 1; sp != SPRAY_FUEL_NUM; ++sp)
            host_real_attribs[PeleLM::pstateY + sp].push_back(0.);
#else
          p.rdata(PeleLM::pstateT) = part_temp;
          p.rdata(PeleLM::pstateDia) = cur_dia;
          p.rdata(PeleLM::pstateRho) = part_rho;
          for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp)
            p.rdata(PeleLM::pstateY + sp) = 0.;
          p.rdata(PeleLM::pstateY) = 1.;
#endif
          host_particles.push_back(p);
          Real pmass = Pi_six*part_rho*std::pow(cur_dia, 3);
          total_mass += num_ppp*pmass;
        }
      }
    }
    auto& particle_tile = GetParticles(lev)[std::make_pair(mfi.index(),
                                                           mfi.LocalTileIndex())];
    auto old_size = particle_tile.GetArrayOfStructs().size();
    auto new_size = old_size + host_particles.size();
    particle_tile.resize(new_size);

    Gpu::copy(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
              particle_tile.GetArrayOfStructs().begin() + old_size);
#ifdef USE_SPRAY_SOA
    for (int i = 0; i != NAR_SPR; ++i) {
      Gpu::copy(Gpu::hostToDevice, host_real_attribs[i].begin(), host_real_attribs[i].end(),
                particle_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
    }
#endif
  }
  // Redistribute is done outside of this function
  return true;
}
  
void
SprayParticleContainer::InitSprayParticles(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
  // Start without any particles
  return;
}
