#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <WritePlotFile.H>
#include <AMReX_AmrData.H>

#include <turbinflow.H>

using namespace amrex;

void PlotFileFromMF(const MultiFab& mf,
                    const Geometry& geom,
                    const std::string& oFile)
{
  Real bgVal = 0;
  IntVect refRatio(D_DECL(2,2,2));
  Vector<std::string> names(mf.nComp());
  for (int i=0; i<mf.nComp(); ++i) {
    names[i] = amrex::Concatenate("MF_",i,5);
  }
    
  writePlotFile(oFile.c_str(),mf,geom,refRatio,bgVal,names);
}

static
void
Extend (FArrayBox& xfab,
        FArrayBox& vfab,
        const Box& domain)
{
  Box tbx = vfab.box();

  tbx.setBig(0, domain.bigEnd(0) + 3);

  const int ygrow = BL_SPACEDIM==3 ? 3 : 1;

  tbx.setBig(1, domain.bigEnd(1) + ygrow);

  xfab.resize(tbx,1);

  Box orig_box = vfab.box();
  vfab.shift(0, 1);
  vfab.shift(1, 1);
  xfab.copy(vfab); // (0,0)

  vfab.shift(0, domain.length(0)); // (1,0)
  xfab.copy(vfab);
  vfab.shift(1, domain.length(1)); // (1,1)
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0)); // (0,1)
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0)); // (-1,1)
  xfab.copy(vfab);
  vfab.shift(1, -domain.length(1)); // (-1,0)
  xfab.copy(vfab);
  vfab.shift(1, -domain.length(1)); // (-1,-1)
  xfab.copy(vfab);
  vfab.shift(0, domain.length(0)); // (0,-1)
  xfab.copy(vfab);
  vfab.shift(0, domain.length(0)); // (1,-1)
  xfab.copy(vfab);
  vfab.shift(0, -domain.length(0) - 1);
  vfab.shift(1,  domain.length(1) - 1);

  if (vfab.box() != orig_box) Abort("Oops, something bad happened");
}

void
WriteTurb(const std::string& turbfile,
          const MultiFab& data,
          const Box& box_turb,
          const RealBox& rb_turb)
{
  if (ParallelDescriptor::IOProcessor())
    if (!UtilCreateDirectory(turbfile, 0755))
      CreateDirectoryFailed(turbfile);
  
  std::string Hdr = turbfile; Hdr += "/HDR";
  std::string Dat = turbfile; Dat += "/DAT";

  std::ofstream ifsd, ifsh;

  if (ParallelDescriptor::IOProcessor()) {
    ifsh.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsh.good())
      FileOpenFailed(Hdr);

    ifsd.open(Dat.c_str(), std::ios::out|std::ios::trunc);
    if (!ifsd.good())
      FileOpenFailed(Dat);

    int coord_turb(0);
    Array<int,BL_SPACEDIM> per_turb = {D_DECL(1,1,1)};
    Geometry geom_turb(box_turb,rb_turb,coord_turb,per_turb);
    const Real* dx_turb = geom_turb.CellSize();

    //
    // Write the first part of the Turb header.
    // Note that this is solely for periodic style inflow files.
    //
    Box box_turb_io(box_turb);
    box_turb_io.setBig(0, box_turb.bigEnd(0) + 3);
    box_turb_io.setBig(1, box_turb.bigEnd(1) + 3);
    box_turb_io.setBig(2, box_turb.bigEnd(2) + 1);

    ifsh << box_turb_io.length(0) << ' '
         << box_turb_io.length(1) << ' '
         << box_turb_io.length(2) << '\n';

    ifsh << rb_turb.length(0) + 2*dx_turb[0] << ' '
         << rb_turb.length(1) + 2*dx_turb[1] << ' '
         << rb_turb.length(2)                << '\n';

    ifsh << per_turb[0] << ' ' << per_turb[1] << ' ' << per_turb[2] << '\n';
  }

  IntVect sm = box_turb.smallEnd();
  IntVect bg = box_turb.bigEnd();
  int dir = BL_SPACEDIM - 1;
  FArrayBox xfab,TMP;
  //
  // We work on one cell wide Z-planes.
  // We first do the lo BL_SPACEDIM plane.
  // And then all the other planes in xhi -> xlo order.
  //
  for (int d = 0; d < BL_SPACEDIM; ++d)
  {
    bg[dir] = sm[dir];
    {
      Box bx(sm,bg);
      TMP.resize(bx,1);
      data.copyTo(TMP,d,0,1);
      Extend(xfab, TMP, box_turb);
      if (ParallelDescriptor::IOProcessor()) {
        ifsh << ifsd.tellp() << std::endl;
        xfab.writeOn(ifsd);
      }
    }
    for (int i = box_turb.bigEnd(dir); i >= box_turb.smallEnd(dir); i--)
    {
      sm[dir] = i;
      bg[dir] = i;
      Box bx(sm,bg);
      TMP.resize(bx,1);
      data.copyTo(TMP,d,0,1);
      Extend(xfab, TMP, box_turb);
      if (ParallelDescriptor::IOProcessor()) {
        ifsh << ifsd.tellp() << std::endl;
        xfab.writeOn(ifsd);
      }
    }
  }
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string pltfile("plt");  pp.query("plotfile",pltfile);
    
    std::string TurbDir("Turb");

    if (pp.countval("TurbDir")>0) {
       pp.query("TurbDir",TurbDir);
    }
    else if (pp.countval("turb_pltfile")) {
      std::string turb_pltfile;
      pp.get("turb_pltfile",turb_pltfile);
      
      DataServices dataServices(turb_pltfile, Amrvis::NEWPLT);
      if( ! dataServices.AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
        // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
      }
      AmrData& amrData = dataServices.AmrDataRef();

      int idV = amrData.StateNumber("x_velocity");
      int finestLevel = amrData.FinestLevel();
      const auto& xvel = amrData.GetGrids(finestLevel,idV);
      MultiFab pltvel(xvel.boxArray(),xvel.DistributionMap(),AMREX_SPACEDIM,0);
      MultiFab::Copy(pltvel,xvel,0,0,1,0);
      for (int i=1; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Copy(pltvel,amrData.GetGrids(finestLevel,idV+i),0,i,1,0);
      }
      
      WriteTurb(TurbDir,pltvel,amrData.ProbDomain()[finestLevel],
                RealBox(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0])));

    } else
    { 
      if (ParallelDescriptor::IOProcessor())
        if (!UtilCreateDirectory(TurbDir, 0755))
          CreateDirectoryFailed(TurbDir);

      std::string Hdr = TurbDir; Hdr += "/HDR";
      std::string Dat = TurbDir; Dat += "/DAT";

      std::ofstream ifsd, ifsh;

      ifsh.open(Hdr.c_str(), std::ios::out|std::ios::trunc);
      if (!ifsh.good())
        FileOpenFailed(Hdr);

      ifsd.open(Dat.c_str(), std::ios::out|std::ios::trunc);
      if (!ifsd.good())
        FileOpenFailed(Dat);

      Box box_turb(IntVect(D_DECL(0,0,0)),
                   IntVect(D_DECL(63,63,63)));
      RealBox rb_turb({D_DECL(0,0,0)},
                      {D_DECL(1,1,1)});
      int coord_turb(0);
      Array<int,BL_SPACEDIM> per_turb = {D_DECL(1,1,1)};
      Geometry geom_turb(box_turb,rb_turb,coord_turb,per_turb);
      const Real* dx_turb = geom_turb.CellSize();

      //
      // Write the first part of the Turb header.
      // Note that this is solely for periodic style inflow files.
      //
      Box box_turb_io(box_turb);
      box_turb_io.setBig(0, box_turb.bigEnd(0) + 3);
      box_turb_io.setBig(1, box_turb.bigEnd(1) + 3);
      box_turb_io.setBig(2, box_turb.bigEnd(2) + 1);

      ifsh << box_turb_io.length(0) << ' '
           << box_turb_io.length(1) << ' '
           << box_turb_io.length(2) << '\n';

      ifsh << rb_turb.length(0) + 2*dx_turb[0] << ' '
           << rb_turb.length(1) + 2*dx_turb[1] << ' '
           << rb_turb.length(2)                << '\n';

      ifsh << per_turb[0] << ' ' << per_turb[1] << ' ' << per_turb[2] << '\n';

      // Create a field to shove in
      FArrayBox vel_turb(box_turb,BL_SPACEDIM);
      Array4<Real> const& fab = vel_turb.array();
      AMREX_PARALLEL_FOR_3D ( box_turb, i, j, k,
                              {
                                Real x = (i+0.5)*dx_turb[0] + rb_turb.lo()[0];
                                Real y = (j+0.5)*dx_turb[1] + rb_turb.lo()[1];
                                Real z = (k+0.5)*dx_turb[2] + rb_turb.lo()[2];
                                Real blobr = 0.25 * (rb_turb.hi()[0] - rb_turb.lo()[0]);
                                Real blobx = 0.5 * (rb_turb.hi()[0] + rb_turb.lo()[0]);
                                Real bloby = 0.5 * (rb_turb.hi()[0] + rb_turb.lo()[0]);
                                Real blobz = 0.5 * (rb_turb.hi()[0] + rb_turb.lo()[0]);
                                Real r = std::sqrt((x-blobx)*(x-blobx) +
                                                   (y-bloby)*(y-bloby) +
                                                   (z-blobz)*(z-blobz));
                                fab(i,j,k,0) = r <= blobr ? 1 : 0;
                                fab(i,j,k,1) = 2.;
                                fab(i,j,k,2) = 3.;
                              });

      // Write "turbulence" as a plotfile
      BoxArray ba(box_turb);
      DistributionMapping dm(ba);
      MultiFab vel(ba,dm,BL_SPACEDIM,0);
      if (dm[0]==ParallelDescriptor::MyProc()) {
        vel[0].copy(vel_turb);
      }
      PlotFileFromMF(vel,geom_turb,Concatenate(pltfile,0));

      // Dump field as a "turbulence file"
      IntVect sm = box_turb.smallEnd();
      IntVect bg = box_turb.bigEnd();
      int dir = BL_SPACEDIM - 1;
      FArrayBox xfab,TMP;
      //
      // We work on one cell wide Z-planes.
      // We first do the lo BL_SPACEDIM plane.
      // And then all the other planes in xhi -> xlo order.
      //
      for (int d = 0; d < BL_SPACEDIM; ++d)
      {
        bg[dir] = sm[dir];
        {
          Box bx(sm,bg);
          TMP.resize(bx,1);
          TMP.copy(vel_turb,bx,d,bx,0,1);
          Extend(xfab, TMP, box_turb);
          ifsh << ifsd.tellp() << std::endl;
          xfab.writeOn(ifsd);
        }
        for (int i = box_turb.bigEnd(dir); i >= box_turb.smallEnd(dir); i--)
        {
          sm[dir] = i;
          bg[dir] = i;
          Box bx(sm,bg);
          TMP.resize(bx,1);
          TMP.copy(vel_turb,bx,d,bx,0,1);
          Extend(xfab, TMP, box_turb);
          ifsh << ifsd.tellp() << std::endl;
          xfab.writeOn(ifsd);
        }
      }
    }

    // Now that turbulence file written, read from it to fill the result
    Real loc_scale = 1;
    Real vel_scale = 1;
    pp.query("loc_scale",loc_scale);
    pp.query("vel_scale",vel_scale);
    TurbParm tp;
    tp.tph = new TurbParmHost();
    init_turbinflow(TurbDir,loc_scale,vel_scale,{{0.,0.}},tp);

    Box box_res(IntVect(D_DECL(0,0,0)),
                IntVect(D_DECL(255,255,127)));
    RealBox rb_res({D_DECL(-1,-1,0)},
                   {D_DECL(1,1,1)});
    pp.query("box_res",box_res);
    Vector<Real> dlo(BL_SPACEDIM), dhi(BL_SPACEDIM);
    if (pp.countval("dlo")>0) {
      pp.getarr("dlo",dlo,0,BL_SPACEDIM);
      pp.getarr("dhi",dhi,0,BL_SPACEDIM);
      rb_res = RealBox(&(dlo[0]),&(dhi[0]));
    }
    int coord_res(0);
    Array<int,BL_SPACEDIM> per_res = {D_DECL(1,1,1)};
    Geometry geom_res(box_res,rb_res,coord_res,per_res);
    Box domain = geom_res.Domain();
    BoxArray ba_res(domain);
    ba_res.maxSize(16);
    DistributionMapping dmap_res(ba_res);
    MultiFab mf_res(ba_res,dmap_res,BL_SPACEDIM,0);

    mf_res.setVal(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf_res,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      fill_with_turb(mfi.tilebox(),mf_res[mfi],0,geom_res,tp);
    }

    std::string outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
    PlotFileFromMF(mf_res,geom_res,outfile);
  }
  Finalize();
}
