#include <AMReX_REAL.H>

module probdata_module

  use mod_Fvar_def, only: maxspec

  implicit none
  
  REAL_T vorterr, temperr, adverr, tempgrad, flametracval, twall
  REAL_T hrr_thresh, RO2_thresh, OH_thresh, mix_thresh, zmax_mix, diff_thresh
  REAL_T splitx, splity, traceSpecVal, xfrontw, yfrontw, zmax_diff
  REAL_T standoff
  REAL_T domnlo(3), domnhi(3)

  integer max_vort_lev, max_temp_lev, max_trac_lev, max_nozzle_lev
  integer max_mix_lev, max_hrr_lev, max_diff_lev
  integer fuelID, oxidID, prodID, nspecies
  character*50 probtype

  integer refine_nozzle
  REAL_T refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, blobx, bloby, blobz, blobr, blobT, Tfrontw, xcen

  REAL_T v_strength,v_width,v_xcen,v_ycen, v_cl_x

  REAL_T v_blob_r, v_blob_T, v_blob_airfrac, turb_scale,phi_in

  REAL_T stTh, Rf, V_in, V_co, T_in, T_co, R_hot, R_hotBL, stBL, pipeTh, pipeBL
  REAL_T tV_in_l, tV_in_r, V_in_new, tV_co_l, tV_co_r, V_co_new, H2_frac
  REAL_T T_stick, rho_stick, h_stick, fuel_N2_vol_percent,T_switch,V_jet
  REAL_T pertmag

  integer :: iN2 = -1
  integer :: iO2 = -1
  integer :: iCH4 = -1
  integer :: iH2 = -1
  integer :: iNC12H26 = -1
  
  integer dimFile
  logical forceInflow
  REAL_T convVelInit, probSizeFile(3), dxFile(3)

  logical bcinit
  integer, parameter :: Nzones=5
  REAL_T u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
  REAL_T Y_bc(0:maxspec-1, Nzones), T_bc(Nzones)
  REAL_T h_bc(Nzones)
  
  integer, parameter :: BL_FUELPIPE = 1
  integer, parameter :: BL_COFLOW   = 2
  integer, parameter :: BL_AMBIENT  = 3
  integer, parameter :: BL_STICK    = 4
  integer, parameter :: BL_WALL     = 5
  integer, parameter :: BL_VOLUME   = 6
  integer, parameter :: BL_PIPEEND  = 7
  
  integer MAXPNTS
  parameter(MAXPNTS = 50)
  REAL_T time_points(0:MAXPNTS),vel_points(0:MAXPNTS),cntl_points(0:MAXPNTS)

  character(50) ac_hist_file
  REAL_T tau_control, cfix, coft_old, sest, V_in_old, corr, &
       changeMax_control, tbase_control, dV_control, scale_control, &
       zbase_control, h_control, controlVelMax
  integer navg_pnts
  integer pseudo_gravity

contains

!subroutines here

end module probdata_module
