#include "prob.H"
#include "turbinflow.H"
#include "AMReX_Utility.H"

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  // amrex::Print() << "myfile="<<myfile<<" \n";

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  // amrex::Print() << "read_pmf2! \n";

  amrex::Vector<std::string> pmf_names;
  pmf_names.resize(variable_count);
  pos1 = 0;
  // pos2 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pmf_names[i] = firstline.substr(pos1 + 1, pos2 - (pos1 + 1));
    pos1 = pos2 + 1;
  }

  // amrex::Print() << "read_pmf3! \n";

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::h_prob_parm_device->pmf_N = line_count;
  PeleC::h_prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    PeleC::h_prob_parm_device->pmf_N * PeleC::h_prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    PeleC::h_prob_parm_device->pmf_N * PeleC::h_prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (unsigned int i = 0; i < PeleC::h_prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (unsigned int j = 0; j < PeleC::h_prob_parm_device->pmf_M; j++) {
      sinput >> PeleC::prob_parm_host
                  ->h_pmf_Y[j * PeleC::h_prob_parm_device->pmf_N + i];
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
  PeleC::d_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::d_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real molefrac[NUM_SPECIES];
  amrex::Real massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};

  if (PeleC::h_prob_parm_device->phi_in < 0) {
    const amrex::Real yl = 0.0;
    const amrex::Real yr = 0.0;
    pmf(yl, yr, pmf_vals, *PeleC::h_prob_parm_device);
    amrex::Real mysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = amrex::max<amrex::Real>(0.0, pmf_vals[3 + n]);
      mysum += molefrac[n];
    }
    molefrac[N2_ID] = 1.0 - (mysum - molefrac[N2_ID]);
    T = pmf_vals[0];
    PeleC::h_prob_parm_device->vn_in = pmf_vals[1];
  } else {
    const amrex::Real a = 0.5;
    for (amrex::Real& n : molefrac) {
      n = 0.0;
    }
    molefrac[O2_ID] =
      1.0 / (1.0 + PeleC::h_prob_parm_device->phi_in / a + 0.79 / 0.21);
    molefrac[H2_ID] = PeleC::h_prob_parm_device->phi_in * molefrac[O2_ID] / a;
    molefrac[N2_ID] = 1.0 - molefrac[H2_ID] - molefrac[O2_ID];
    T = PeleC::h_prob_parm_device->T_in;
  }
  const amrex::Real p = PeleC::h_prob_parm_device->pamb;

  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::h_prob_parm_device->fuel_state[URHO] = rho;
  PeleC::h_prob_parm_device->fuel_state[UMX] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UMY] = rho * vt;
  PeleC::h_prob_parm_device->fuel_state[UMZ] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UEINT] = rho * e;
  PeleC::h_prob_parm_device->fuel_state[UEDEN] = rho * (e + ek);
  PeleC::h_prob_parm_device->fuel_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::h_prob_parm_device->fuel_state[UFS + n - 1] = rho * massfrac[n];
  }

  //  amrex::Print() << " rho=" << rho  << " \n";
  // amrex::Print() << " vt =" << vt  << " \n";
  // amrex::Print() << " e  =" << e  << " \n";
  // amrex::Print() << " ek =" << ek  << " \n";
  // amrex::Print() << " T  =" << T  << " \n";
  // for (int n = 0; n < NUM_SPECIES; n++) {
  //  amrex::Print() << " Y" << "("<<n<<") =" <<  massfrac[n]   << " \n";
  //}
}

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  std::string pmf_datafile;

  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);
  pp.query("phi_in", PeleC::h_prob_parm_device->phi_in);
  pp.query("pertmag", PeleC::h_prob_parm_device->pertmag);
  pp.query("T_inlet_pert_mag", PeleC::h_prob_parm_device->T_inlet_pert_mag);
  pp.query("T_inlet_pert_freq", PeleC::h_prob_parm_device->T_inlet_pert_freq);
  pp.query("P_inlet_pert_mag", PeleC::h_prob_parm_device->P_inlet_pert_mag);
  pp.query("P_inlet_pert_freq", PeleC::h_prob_parm_device->P_inlet_pert_freq);
  pp.query("pmf_datafile", pmf_datafile);

  PeleC::h_prob_parm_device->L[0] = probhi[0] - problo[0];
  PeleC::h_prob_parm_device->L[1] = probhi[1] - problo[1];
  PeleC::h_prob_parm_device->L[2] = probhi[2] - problo[2];

  read_pmf(pmf_datafile);

  PeleC::prob_parm_host->do_turb = false;
  if (pp.countval("turb_file") > 0) {
#if AMREX_SPACEDIM==2
    amrex::Abort("Turbulence inflow unsupported in 2D.");
#endif
  std::string turb_file = "";
    pp.query("turb_file", turb_file);
    amrex::Real turb_scale_loc = 1.0;
    pp.query("turb_scale_loc", turb_scale_loc);
    amrex::Real turb_scale_vel = 1.0;
    pp.query("turb_scale_vel", turb_scale_vel);

    PeleC::prob_parm_host->do_turb = true;

    // Hold nose here - required because of dynamically allocated data in tp
    AMREX_ASSERT_WITH_MESSAGE(
      PeleC::h_prob_parm_device->tp.tph == nullptr,
      "Can only be one TurbParmHost");
    PeleC::h_prob_parm_device->tp.tph = new TurbParmHost;

    amrex::Vector<amrex::Real> turb_center = {
      {0.5 * (probhi[0] + problo[0]), 0.5 * (probhi[1] + problo[1])}};
    pp.queryarr("turb_center", turb_center);
    AMREX_ASSERT_WITH_MESSAGE(
      turb_center.size() == 2, "turb_center must have two elements");
    for (int n = 0; n < turb_center.size(); ++n) {
      turb_center[n] *= turb_scale_loc;
    }
    int turb_nplane = 32;
    pp.query("turb_nplane", turb_nplane);
    AMREX_ASSERT(turb_nplane > 0);
    amrex::Real turb_conv_vel = 1;
    pp.query("turb_conv_vel", turb_conv_vel);
    AMREX_ASSERT(turb_conv_vel > 0);
    init_turbinflow(
      turb_file, turb_scale_loc, turb_scale_vel, turb_center, turb_conv_vel,
      turb_nplane, PeleC::h_prob_parm_device->tp);
  }

  init_bc();
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}

#if AMREX_SPACEDIM==2
#define D_DROP_FRONT(a,b,c) b,c
#elif AMREX_SPACEDIM==3
#define D_DROP_FRONT(a,b,c) a,b,c
#endif

void
EBUTurn::build(const amrex::Geometry& geom, const int max_coarsening_level)
{
  // Under arbeid. Fiks parametrisering og input.
  amrex::ParmParse pp("exp_chan");
  amrex::Vector<amrex::Real> box1lo, box1hi, box2lo, box2hi;
  amrex::Real p1_y, p2_y;
  pp.getarr("b1_corner_lo", box1lo);
  pp.getarr("b1_corner_hi", box1hi);
  pp.getarr("b2_corner_lo", box2lo);
  pp.getarr("b2_corner_hi", box2hi);
  pp.get("exp_y_lo", p1_y);
  pp.get("exp_y_hi", p2_y);

  // p2_y = 1.8;
  // p1_y = 0.0;
  const amrex::Real cen_y = 0.5 * (p2_y - p1_y);
  amrex::Real dx = geom.CellSize(0); // remember dx = dy = dz

  amrex::Print() << "Hellø! \n";
  amrex::Print() << "Wi nøt trei a høliday in Sweden this yër?\n";

  // Main part: Box with half cylinder added to the left
  // Put a margin of 0.3*dx to make sure the EB algorithm finds the boundaries
  // in y and not the z boundary. (Needed?)
  amrex::EB2::BoxIF box1(
    {D_DROP_FRONT(-1.0, 0.0 + 0.3 * dx, cen_y)}, {D_DROP_FRONT(20.0, p2_y, 8.0 + 0.3 * dx)}, true);
  //    EB2::CylinderIF cyl(Real a_radius, int a_direction,
  //            const RealArray& a_center, bool a_inside)
  //    direction == 0 for x etc.
  const int cyl_axis = AMREX_D_PICK(-1, 2, 0);
  amrex::EB2::CylinderIF cyl(cen_y - 0.3 * dx, cyl_axis, {D_DROP_FRONT(0.0, cen_y, cen_y)}, true);
  auto channel = amrex::EB2::makeIntersection(cyl, box1);

  // The upper 'extension' of the channel
  // Straight part
  amrex::EB2::BoxIF box2(
    {D_DROP_FRONT(-1.0, p2_y, 4.2)}, {D_DROP_FRONT(20.0, 2.0 - 0.3 * dx, 8.0 + 0.3 * dx)}, true);
  auto channel2 = amrex::EB2::makeIntersection(channel, box2);
  // The tip.  Construct triangle from planes:
  amrex::EB2::PlaneIF plane1({D_DROP_FRONT(0.0, p2_y, 4.2)}, {D_DROP_FRONT(0.0, -1.0, 0.0)});
  amrex::EB2::PlaneIF plane2({D_DROP_FRONT(0.0, p2_y, 4.2)}, {D_DROP_FRONT(0.0, 0.0, 1.0)});
  amrex::EB2::PlaneIF plane3({D_DROP_FRONT(0.0, p2_y, 3.8)}, {D_DROP_FRONT(0.0, 2.0, -1.0)});
  auto triangle = amrex::EB2::makeUnion(plane1, plane2, plane3);
  // I'm not allowed to say: triangle = amrex::EB2::makeComplement(triangle);
  // auto triangle4 = amrex::EB2::makeComplement(triangle3);
  auto channel3 = amrex::EB2::makeIntersection(channel2, triangle);

  // The separation wall
  // Straight part
  amrex::EB2::BoxIF box3(
    {D_DROP_FRONT(-1.0, cen_y, 1.2)}, {D_DROP_FRONT(20.0, 1.0, 8.0 + 0.3 * dx)}, false);
  auto channel4 = amrex::EB2::makeUnion(channel3, box3);
  // The tip. Construct triangle from planes:
  amrex::EB2::PlaneIF plane4({D_DROP_FRONT(0.0, cen_y, 1.2)}, {D_DROP_FRONT(0.0, -1.0, 0.0)});
  amrex::EB2::PlaneIF plane5({D_DROP_FRONT(0.0, cen_y, 1.2)}, {D_DROP_FRONT(0.0, 0.0, 1.0)});
  amrex::EB2::PlaneIF plane6({D_DROP_FRONT(0.0, cen_y, 1.0)}, {D_DROP_FRONT(0.0, 2.0, -1.0)});
  amrex::EB2::PlaneIF plane7({D_DROP_FRONT(0.0, cen_y, 1.023)}, {D_DROP_FRONT(0.0, 0.0, -1.0)});
  auto triangle3 = amrex::EB2::makeUnion(plane4, plane5, plane6, plane7);
  auto triangle4 = amrex::EB2::makeComplement(triangle3);
  auto channel5 = amrex::EB2::makeUnion(channel4, triangle4);

  auto gshop = amrex::EB2::makeShop(channel5);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}

void
EBUTurnFlipped::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  // In progress. Parameterization and input to be fixed.
  amrex::ParmParse pp("exp_chan");
  amrex::Real p1_y, p2_y;
  pp.get("exp_y_lo", p1_y);
  pp.get("exp_y_hi", p2_y);

  const amrex::Real cen_y = 0.5 * (p2_y - p1_y);
  amrex::Real dx = geom.CellSize(0); // remember dx = dy = dz

  // Main part: Box with half cylinder added to the right
  // Put a margin of 0.3*dx to make sure the EB algorithm finds the boundaries
  // in y and not the z boundary. (Needed?)
  amrex::EB2::BoxIF box1(
    {D_DROP_FRONT(-1.0, 0.0 + 0.3 * dx, -0.3 * dx)}, {D_DROP_FRONT(20.0, p2_y, 8.0 - cen_y)}, true);
  //    EB2::CylinderIF cyl(Real a_radius, int a_direction,
  //            const RealArray& a_center, bool a_inside)
  //    direction == 0 for x osv.
  const int cyl_axis = AMREX_D_PICK(-1, 2, 0);
  amrex::EB2::CylinderIF cyl(
    cen_y - 0.3 * dx, cyl_axis, {D_DROP_FRONT(0.0, cen_y, 8.0 - cen_y)}, true);
  auto channel = amrex::EB2::makeIntersection(cyl, box1);

  // The upper 'extension' of the channel
  // Straight part
  amrex::EB2::BoxIF box2(
    {D_DROP_FRONT(-1.0, p2_y, -0.3 * dx)}, {D_DROP_FRONT(20.0, 2.0 - 0.3 * dx, 3.8)}, true);
  auto channel2 = amrex::EB2::makeIntersection(channel, box2);
  // The tip.  Construct triangle from planes:
  amrex::EB2::PlaneIF plane1({D_DROP_FRONT(0.0, p2_y, 3.8)}, {D_DROP_FRONT(0.0, -1.0, 0.0)});
  amrex::EB2::PlaneIF plane2({D_DROP_FRONT(0.0, p2_y, 3.8)}, {D_DROP_FRONT(0.0, 0.0, -1.0)});
  amrex::EB2::PlaneIF plane3({D_DROP_FRONT(0.0, p2_y, 4.2)}, {D_DROP_FRONT(0.0, 2.0, 1.0)});
  auto triangle = amrex::EB2::makeUnion(plane1, plane2, plane3);
  // I'm not allowed to say: triangle = amrex::EB2::makeComplement(triangle);
  // auto triangle4 = amrex::EB2::makeComplement(triangle3);
  auto channel3 = amrex::EB2::makeIntersection(channel2, triangle);

  // The separation wall
  // Straight part
  //    amrex::EB2::BoxIF box3({-1.0,cen_y,-0.3*dx}, {20.0,1.0,6.8}, false);
  // make an y coordinate equal to 1.0 - 13*dx
  // for dx = 2 / 256
  const amrex::Real y_sep = 0.8984375;
  amrex::EB2::BoxIF box3({D_DROP_FRONT(-1.0, y_sep, -0.3 * dx)}, {D_DROP_FRONT(20.0, 1.0, 6.8)}, false);
  auto channel4 = amrex::EB2::makeUnion(channel3, box3);
  // The tip. Construct triangle from planes:
  amrex::EB2::PlaneIF plane4({D_DROP_FRONT(0.0, y_sep, 6.8)}, {D_DROP_FRONT(0.0, -1.0, 0.0)});
  amrex::EB2::PlaneIF plane5({D_DROP_FRONT(0.0, y_sep, 6.8)}, {D_DROP_FRONT(0.0, 0.0, -1.0)});
  amrex::EB2::PlaneIF plane6({D_DROP_FRONT(0.0, y_sep, 7.0)}, {D_DROP_FRONT(0.0, 2.0, 1.0)});
  amrex::EB2::PlaneIF plane7({D_DROP_FRONT(0.0, y_sep, 6.977)}, {D_DROP_FRONT(0.0, 0.0, 1.0)});
  auto triangle3 = amrex::EB2::makeUnion(plane4, plane5, plane6, plane7);
  auto triangle4 = amrex::EB2::makeComplement(triangle3);
  auto channel5 = amrex::EB2::makeUnion(channel4, triangle4);

  auto gshop = amrex::EB2::makeShop(channel5);
  amrex::EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
}
