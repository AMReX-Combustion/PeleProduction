#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

#include "PeleC.H"
#include "prob.H"
#include "turbinflow.H"
#include "PelePhysics.H"
#include "IndexDefines.H"

struct PCHypFillExtDir
{
  ProbParmDevice const* lprobparm;

  AMREX_GPU_HOST
  constexpr explicit PCHypFillExtDir(const ProbParmDevice* d_prob_parm)
    : lprobparm(d_prob_parm)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    // const amrex::Real* prob_hi = geom.ProbHi();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

    const int* bc = bcr->data();

    amrex::Real s_int[NVAR] = {0.0};
    amrex::Real s_ext[NVAR] = {0.0};

    // xlo and xhi
    int idir = 0;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domlo[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(domhi[idir], iv[1], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#if AMREX_SPACEDIM > 1
    // ylo and yhi
    idir = 1;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domlo[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      amrex::IntVect loc(AMREX_D_DECL(iv[0], domhi[idir], iv[2]));
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(loc, n);
      }
      bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#if AMREX_SPACEDIM == 3
    // zlo and zhi
    idir = 2;
    if ((bc[idir] == amrex::BCType::ext_dir) && (iv[idir] < domlo[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domlo[idir], n);
      }
      bcnormal(x, s_int, s_ext, idir, +1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    } else if (
      (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) &&
      (iv[idir] > domhi[idir])) {
      for (int n = 0; n < NVAR; n++) {
        s_int[n] = dest(iv[0], iv[1], domhi[idir], n);
      }
      bcnormal(x, s_int, s_ext, idir, -1, time, geom, *lprobparm);
      for (int n = 0; n < NVAR; n++) {
        dest(iv, n) = s_ext[n];
      }
    }
#endif
#endif
  }
};

struct PCReactFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dest*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
  }
};

void
StateToPrim(const amrex::FArrayBox& state,
            const amrex::Box&       box,
            amrex::FArrayBox&       prim,
            int                     clean_massfrac)
{
  const auto& astate = state.array();
  const auto& aprim = prim.array();
  amrex::ParallelFor(box, [astate, aprim, clean_massfrac]
  AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real rho = astate(i,j,k,URHO);
    amrex::Real rhoinv = 1.0 / rho;
    aprim(i,j,k,QRHO) = rho;
    aprim(i,j,k,QU) = astate(i,j,k,UMX) * rhoinv;
    aprim(i,j,k,QV) = astate(i,j,k,UMY) * rhoinv;
    aprim(i,j,k,QW) = astate(i,j,k,UMZ) * rhoinv;
    amrex::Real kineng = 0.5 * rho * (aprim(i,j,k,QU)*aprim(i,j,k,QU)
                                      + aprim(i,j,k,QV)*aprim(i,j,k,QV)
                                      + aprim(i,j,k,QW)*aprim(i,j,k,QW));
    aprim(i,j,k,QREINT) = astate(i,j,k,UEDEN) - kineng;
  });
}

void
PrimToState(const amrex::FArrayBox& prim,
            const amrex::Box&       box,
            amrex::FArrayBox&       state,
            int                     clean_massfrac)
{
  const auto& aprim = prim.array();
  const auto& astate = state.array();
  amrex::ParallelFor(box, [aprim, astate, clean_massfrac]
  AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real rho = aprim(i,j,k,URHO);
    amrex::Real rhoinv = 1.0 / rho;
    astate(i,j,k,QRHO) = rho;
    astate(i,j,k,UMX) = aprim(i,j,k,QU) * rho;
    astate(i,j,k,UMY) = aprim(i,j,k,QV) * rho;
    astate(i,j,k,UMZ) = aprim(i,j,k,QW) * rho;
    amrex::Real kineng = 0.5 * rho * (aprim(i,j,k,QU)*aprim(i,j,k,QU)
                                      + aprim(i,j,k,QV)*aprim(i,j,k,QV)
                                      + aprim(i,j,k,QW)*aprim(i,j,k,QW));

    astate(i,j,k,UEDEN) = aprim(i,j,k,QREINT) + kineng;
  });
}

void
pc_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;
  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(
    PCHypFillExtDir{lprobparm});
  hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);

  if (PeleC::prob_parm_host->do_turb) {

    int clean_massfrac = 0; // FIXME: get this from PeleC

    AMREX_ASSERT_WITH_MESSAGE(scomp==0 && numcomp==NVAR,"Fluctations code requires group bc filler approach");
    for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
      auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(),dir) & bx);
      auto isect_lo = bndryBoxLO.ok();
      if (bcr[1].lo()[dir]==EXT_DIR && isect_lo) {
        amrex::FArrayBox vel_fluctsLO(bndryBoxLO,QVAR);
        StateToPrim(data,bndryBoxLO,vel_fluctsLO,clean_massfrac);
        add_turb(bx, vel_fluctsLO, QU, geom, time, dir, amrex::Orientation::low, PeleC::d_prob_parm_device->tp);
        PrimToState(vel_fluctsLO,bndryBoxLO,data,clean_massfrac);
      }

      auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(),dir) & bx);
      auto isect_hi = bndryBoxHI.ok();
      if (bcr[1].hi()[dir]==EXT_DIR && isect_hi) {
        amrex::FArrayBox vel_fluctsHI(bndryBoxHI,QVAR);
        StateToPrim(data,bndryBoxHI,vel_fluctsHI,clean_massfrac);
        add_turb(bx, vel_fluctsHI, QU, geom, time, dir, amrex::Orientation::high, PeleC::d_prob_parm_device->tp);
        PrimToState(vel_fluctsHI,bndryBoxHI,data,clean_massfrac);
      }
    }
  }
}

#ifdef PELEC_USE_REACTIONS
void
pc_reactfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  amrex::GpuBndryFuncFab<PCReactFillExtDir> react_bndry_func(
    PCReactFillExtDir{});
  react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
#endif

void
pc_nullfill(
  amrex::Box const& /*bx*/,
  amrex::FArrayBox& /*data*/,
  const int /*dcomp*/,
  const int /*numcomp*/,
  amrex::Geometry const& /*geom*/,
  const amrex::Real /*time*/,
  const amrex::Vector<amrex::BCRec>& /*bcr*/,
  const int /*bcomp*/,
  const int /*scomp*/)
{
}
