#include <AMReX_REAL.H>
#include "mechanism.h"

module probdata_module

  use amrex_error_module, only : amrex_abort
  implicit none

  integer, parameter :: nspecies = NUM_SPECIES

  REAL_T vorterr, temperr, adverr, tempgrad, flametracval, twall
  REAL_T jet_mixfrac, jet_hrr, jet_zmax, &
         jet_core_mixfrac, jet_core_zmax, &
         ltc_RO2, &
         edge_flame_OH, edge_flame_zmax, &
         prmx_flame_HO2, prmx_flame_mixfrac, prmx_flame_zmax
  REAL_T splitx, splity, traceSpecVal, xfrontw, yfrontw, zmax_diff, zmax_edge
  REAL_T standoff
  REAL_T domnlo(3), domnhi(3)

  integer max_vort_lev, max_temp_lev, max_trac_lev, max_nozzle_lev
  integer max_mix_lev, max_hrr_lev, max_diff_lev
  integer fuelID, oxidID, prodID

  integer refine_nozzle
  REAL_T refine_nozzle_x, refine_nozzle_y, refine_nozzle_z, blobx, bloby, blobz, blobr, blobT, Tfrontw, xcen

  REAL_T v_strength,v_width,v_xcen,v_ycen, v_cl_x

  REAL_T v_blob_r, v_blob_T, v_blob_airfrac
  REAL_T T_ox, T_fu, mixfrac_in, mixfrac_co
  logical h_fu_nist

  REAL_T stTh, Rf, V_in, V_co, R_hot, R_hotBL, stBL, pipeTh, pipeBL
  REAL_T tV_in_l, tV_in_r, V_in_new, tV_co_l, tV_co_r, V_co_new, H2_frac
  REAL_T T_stick, rho_stick, h_stick, T_switch
  REAL_T pertmag

  integer :: iN2 = -1

  logical :: do_flct = .false.
  REAL_T  :: turb_scale = 1.d0
  REAL_T probSizeFile(3), dxFile(3)

  REAL_T Y_fu(nspecies), Y_ox(nspecies)

  logical bcinit
  integer, parameter :: Nzones=5
  REAL_T u_bc(Nzones), v_bc(Nzones), w_bc(Nzones), rho_bc(Nzones)
  REAL_T Y_bc(0:nspecies-1, Nzones), T_bc(Nzones)
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

  public :: get_nist_enthalpy, get_T_from_hY


  contains

    subroutine get_T_from_hY(T, h, Y)
      use PeleLM_F,  only: pphys_TfromHYpt

      implicit none
      REAL_T,                      intent(inout) :: T
      REAL_T,                      intent(in) :: h
      REAL_T, dimension(nspecies), intent(in) :: Y

      REAL_T,  parameter :: max_err = BL_REAL_E(3.9,-11)
      integer, parameter :: max_iter = 200
      REAL_T res(0:max_iter-1)
      integer N

      call pphys_TfromHYpt(T, h, Y, max_err, max_iter, res, N)

      if (N .lt. 0) then
        write(6,'(a)')    'pphys_TfromHYpt:'
        write(6,'(a,f7.2)')  '  T        = ', T
        write(6,'(a,f9.2)')  '  h        = ', h
        write(6,*)           '  Y        = ', Y
        write(6,'(a,i3)')    '  max_iter = ', max_iter
        write(6,'(a,e15.7)') '  max_err  = ', max_err
        write(6,'(a,i3)')    '  Niter    = ', N
        write(6,*)           '  res      = ', res
        call amrex_abort('Something went wrong in get_T_from_hY')
      end if
    end subroutine get_T_from_hY


    function get_nist_enthalpy(name, temp) result(h)
      !-------------------------------------------------------------------------
      !
      ! Calculate the enthalpy of the fuel based on the liquid properties
      ! obtained from NIST. This thus accounts for the heat of evaporation.
      ! Currently only implemented for n-dodecane, n-heptane and iso-octane.
      !
      ! IN:
      !   name: fuel name
      !   temp: temperature (K)
      !
      ! OUT:
      !   h: enthalpy (J/(kg.K))
      !
      ! 1.  Liquid n-dodecane data from NIST:
      !     https://webbook.nist.gov/cgi/cbook.cgi?ID=C112403&Mask=2#Thermo-Condensed
      !
      !     W                      (kg/kmol)   = 170.3348
      !     DeltaH_liq             (kJ/mol)    = -352.1 +/- 1.4
      !     C_{p,liq} @ T=298.15 K (J/(mol.K)) = 376.00
      !
      ! 2.  Liquid n-heptane data from NIST:
      !     https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=2#Thermo-Condensed
      !
      !     W                      (kg/kmol)   = 100.2019
      !     DeltaH_liq             (kJ/mol)    = -224.4 +/- 0.79
      !     C_{p,liq} @ T=298.15 K (J/(mol.K)) = 224.64
      !
      ! 3.  Liquid iso-octane data from NIST:
      !     https://webbook.nist.gov/cgi/cbook.cgi?ID=C540841&Mask=2#Thermo-Condensed
      !
      !     W                      (kg/kmol)   = 114.2285
      !     DeltaH_liq             (kJ/mol)    = -259.3 +/- 1.3
      !     C_{p,liq} @ T=298.15 K (J/(mol.K)) = 242.49
      !-------------------------------------------------------------------------

      implicit none
      character*16, intent(in) :: name
      REAL_T, intent(in) :: temp

      REAL_T, parameter :: temp_nist_ref = 298.15 ! K
      REAL_T :: h, W, dh_liq, cp_liq

      if ((index(name,'c12h26').gt.0) .or. (index(name,'C12H26').gt.0)) then
        W          = 170.3348e-3 ! kg/mol
        dh_liq     = -352.1e3    ! J/mol
        cp_liq     = 376.00      ! J/((mol.K))
      else if ((index(name,'c7h16').gt.0) .or. (index(name,'C7H16').gt.0)) then
        W          = 100.2019e-3 ! kg/mol
        dh_liq     = -224.4e3    ! J/mol
        cp_liq     = 224.64      ! J/(mol.K)
      else if ((index(name,'c8h18').gt.0) .or. (index(name,'C8H18').gt.0)) then
        W          = 114.2285e-3 ! kg/mol
        dh_liq     = -259.3e3    ! J/mol
        cp_liq     = 242.49      ! J/(mol.K)
      else
        call amrex_abort('unknown fuel in get_nist_enthalpy')
      endif

      dh_liq = dh_liq/W ! J/kg
      cp_liq = cp_liq/W ! J/(kg.K)
      h = dh_liq + cp_liq*(temp - temp_nist_ref) ! J/kg

      return
    end function get_nist_enthalpy

end module probdata_module
