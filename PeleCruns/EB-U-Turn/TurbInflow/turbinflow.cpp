#include <turbinflow.H>

void
init_turbinflow(const std::string& turb_file,
                amrex::Real        turb_scale_loc,
                amrex::Real        turb_scale_vel,
                TurbParm&          tp)
{
  amrex::Print() << "Initializing turbulence file: " << turb_file
                 << " (location coordinates in will be scaled by " << turb_scale_loc
                 << " and velocity out to be scaled by " << turb_scale_vel << ")" << std::endl;

  tp.tph->turb_file = turb_file;
  tp.nplane = 32;
  tp.turb_scale_loc = turb_scale_loc;
  tp.turb_scale_vel = turb_scale_vel;

  std::string turb_header = tp.tph->turb_file + "/HDR";
  std::ifstream is(turb_header.c_str());
  amrex::Array<int, AMREX_SPACEDIM> npts = {{0}};
  amrex::Array<amrex::Real, AMREX_SPACEDIM> probsize = {{0}};
  amrex::Array<int, AMREX_SPACEDIM> iper = {{0}};
  is >> npts[0] >> npts[1] >> npts[2];
  is >> probsize[0] >> probsize[1] >> probsize[2];
  is >> iper[0] >> iper[1] >> iper[2]; // Unused - we assume it is always fully periodic

  for (int n=0; n<AMREX_SPACEDIM; ++n) {
    tp.dx[n] = probsize[n] / amrex::Real(npts[n]-1);
    tp.dxinv[n] = 1.0/tp.dx[n];
  }

  // one ghost point on each side, tangential to inflow face
  tp.pboxsize[0] = probsize[0] - 2.0*tp.dx[0];
  tp.pboxsize[1] = probsize[1] - 2.0*tp.dx[1];
  tp.pboxsize[2] = probsize[2];
  
  tp.npboxcells[0] = npts[0] - 3;
  tp.npboxcells[1] = npts[1] - 3;
  tp.npboxcells[2] = npts[2];

  tp.pboxlo[0] = -0.5 * tp.pboxsize[0];
  tp.pboxlo[1] = -0.5 * tp.pboxsize[1];
  tp.pboxlo[2] = 0.;

  amrex::Box sbx(amrex::IntVect(AMREX_D_DECL(1,1,1)),
                 amrex::IntVect(AMREX_D_DECL(npts[0], npts[1], tp.nplane)));

  tp.sdata = new amrex::FArrayBox(sbx,3);

  tp.kmax = npts[2];

  amrex::Real rdummy;
  if (tp.isswirltype) {
    for (int i = 0; i < tp.kmax; i++) {
      is >> rdummy; // Time for each plane - unused at the moment
    }
  }
  tp.tph->offset_dv.resize(tp.kmax*AMREX_SPACEDIM);
  tp.offset = tp.tph->offset_dv.data();
  tp.offset_size = tp.tph->offset_dv.size();
  for (int i = 0; i < tp.offset_size; i++) {
    is >> tp.offset[i];
  }  
  is.close();

  tp.turbinflow_initialized = true;
}

void
read_one_turb_plane(int       iplane,
                    int       k,
                    TurbParm& tp)
{
  //
  // There are AMREX_SPACEDIM * kmax planes of FABs.
  // The first component are in the first kmax planes,
  // the second component in the next kmax planes, ....
  // Note also that both (*plane) and (*ncomp) start from
  // 1 not 0 since they're passed from Fortran.
  //
  std::string turb_data = tp.tph->turb_file + "/DAT";
  std::ifstream ifs(turb_data.c_str());
  
  amrex::Box dstBox = tp.sdata->box();
  dstBox.setSmall(2,iplane);
  dstBox.setBig(  2,iplane);

  amrex::Print() << "Loading plane " << k << " into slot " << iplane << std::endl;
    
  for (int n=0; n<AMREX_SPACEDIM; ++n) {

    const long offset_idx = (iplane - 1) + (n * tp.kmax);
    AMREX_ASSERT_WITH_MESSAGE(offset_idx < tp.offset_size, "Bad turb fab offset idx");
    
    const long start = tp.offset[offset_idx];

    ifs.seekg(start, std::ios::beg);

    if (!ifs.good())
      amrex::Abort("getplane(): seekg() failed");
  
    amrex::FArrayBox tmp;
    tmp.readFrom(ifs);
    amrex::Box srcBox = tmp.box();

    tp.sdata->copy(tmp,srcBox,0,dstBox,n,1);
  }
  ifs.close();
}

void
read_turb_planes(amrex::Real z,
                 TurbParm&   tp)
{
  int izlo = z * tp.dxinv[2] - 1;
  tp.szlo = izlo * tp.dx[2];
  tp.szhi = tp.szlo
    + (tp.nplane - 1) * tp.dx[2];

  for (int iplane=1; iplane<=tp.nplane; ++iplane) {
    int k = (izlo+iplane-1 % tp.npboxcells[2]) + 1;
    read_one_turb_plane(iplane,k,tp);
  }
  tp.turbinflow_planes_initialized = true;
}

void
fill_turb_plane(const amrex::Vector<amrex::Real>& x,
                const amrex::Vector<amrex::Real>& y,
                amrex::Real                       z,
                amrex::FArrayBox&                 v,
                TurbParm&                         tp)
{
  if ((!tp.turbinflow_planes_initialized) || 
      (z < tp.szlo + 0.5 * tp.dx[2])     ||
      (z > tp.szhi - 0.5 * tp.dx[2]) )
  {
    read_turb_planes(z,tp);
  }

  v.setVal(0);
  return;

  const auto& bx = v.box();
  const auto& vd = v.array();
  const auto& sd = tp.sdata->array();
  amrex::Array<amrex::Real,3> cz;

  amrex::Real zz = (z - tp.szlo) * tp.dxinv[2];
  int k0 = std::round(zz) - 1;
  zz -= amrex::Real(k0);
  cz[0] = 0.5 * (zz-1.0) * (zz - 2.0);
  cz[1] = zz * (2.0 - zz);
  cz[2] = 0.5 * zz * (zz - 1.0);
  k0 += 1; // because the plane numbering is 1-based
  k0 = amrex::min(amrex::max(k0,1),tp.nplane-2);

  amrex::ParallelFor(bx, [x, y, sd, vd, tp, zz, k0, cz]
  AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    amrex::Array<amrex::Real,3> cx, cy, ydata;
    amrex::Array<amrex::Array<amrex::Real,3>,3> zdata;

    for (int n=0; n<3; ++n) {
      amrex::Real xx = (x[i] - tp.pboxlo[0]) * tp.dxinv[0];
      amrex::Real yy = (y[j] - tp.pboxlo[1]) * tp.dxinv[1];
      int i0 = std::round(xx) - 1;
      int j0 = std::round(yy) - 1;
      xx -= amrex::Real(i0);
      yy -= amrex::Real(j0);
      cx[0] = 0.5 * (xx - 1.0) * (xx - 2.0);
      cy[0] = 0.5 * (yy - 1.0) * (yy - 2.0);
      cx[1] = xx * (2.0 - xx);
      cy[1] = yy * (2.0 - yy);
      cx[2] = 0.5 * xx * (xx - 1.0);
      cy[2] = 0.5 * yy * (yy - 1.0);

      i0 = (i0 % tp.npboxcells[0]) + 2; // ! +2 as j0
      j0 = (j0 % 2) + 2;                //+2 because sdataindex starts with 1 and there is a ghost point <---------- FIXME

      for (int ii=0; ii<=2; ++ii) {
        for (int jj=0; jj<=2; ++jj) {
          zdata[ii][jj] = cz[0]*sd(i0+ii,j0+jj,k0  ,n) 
            +             cz[1]*sd(i0+ii,j0+jj,k0+1,n)
            +             cz[2]*sd(i0+ii,j0+jj,k0+2,n);
        }
      }
      for (int ii=0; ii<=2; ++ii) {
        ydata[ii] = cy[0]*zdata[ii][0] + cy[1]*zdata[ii][1] + cy[2]*zdata[ii][2];
      }
      vd(i,j,n) = cx[0]*ydata[0] + cx[1]*ydata[1] + cx[2]*ydata[2];
    }
  });
}

void
add_turb(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp,
  const int dir,
  const amrex::Orientation::Side& side,
  TurbParm& tp)
{
  AMREX_ASSERT_WITH_MESSAGE(dir == 2, "Sadly, the fluctuation code currently only works in the third dimension");
  AMREX_ASSERT(tp.turbinflow_initialized);
  
  amrex::Box bvalsBox = bx;
  int planeLoc = ( side == amrex::Orientation::low
                   ? geom.Domain().smallEnd()[dir]-1
                   : geom.Domain().bigEnd()[dir]  +1 );
  bvalsBox.setSmall(dir,planeLoc);
  bvalsBox.setBig(  dir,planeLoc);

  amrex::FArrayBox v(bvalsBox,3);
  v.setVal(0);

  amrex::Vector<amrex::Real> x(bvalsBox.size()[0]), y(bvalsBox.size()[1]);
  for (int i=bvalsBox.smallEnd()[0]; i<=bvalsBox.bigEnd()[0]; ++i) {
    x[i-bvalsBox.smallEnd()[0]] = (geom.ProbLo()[0] + (i+0.5)*geom.CellSize(0)) * tp.turb_scale_loc;
  }
  for (int j=bvalsBox.smallEnd()[1]; j<=bvalsBox.bigEnd()[1]; ++j) {
    y[j-bvalsBox.smallEnd()[1]] = (geom.ProbLo()[1] + (j+0.5)*geom.CellSize(1)) * tp.turb_scale_loc;
  }

  amrex::Real z = time / tp.turb_conv_vel;
  fill_turb_plane(x, y, z, v, tp);
  v.mult(tp.turb_scale_vel);
}

