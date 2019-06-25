#include <AMReX_REAL.H>

module probdata_module

  use mod_Fvar_def, only: maxspec

  implicit none

  integer, parameter :: BL_FUELPIPE = 1
  integer, parameter :: BL_OUTFLOW  = 2
  integer, parameter :: BL_OXIDIZER = 3
  integer, parameter :: BL_AIR      = 4
  integer, parameter :: BL_PIPEEND  = 5
  integer, parameter :: BL_VOLUME   = 6
  integer, parameter :: BL_NZONES   = 6

  integer, parameter :: BL_XLO = 0
  integer, parameter :: BL_YLO = 1
  integer, parameter :: BL_XHI = 2
  integer, parameter :: BL_YHI = 3
  integer, parameter :: BL_ZLO = 4
  integer, parameter :: BL_ZHI = 5
  integer, parameter :: BL_INTERIOR = 6

  logical :: bcinit = .false.
  
  REAL_T :: u_bc(BL_NZONES), v_bc(BL_NZONES), w_bc(BL_NZONES)
  REAL_T :: Y_bc(0:maxspec-1,BL_NZONES), T_bc(BL_NZONES), h_bc(BL_NZONES), rho_bc(BL_NZONES)

  REAL_T fuel_ox_split, ox_air_split, blobw, blobx, bloby, blobr, blobT, pipeBL, pipeTh,&
       T_air, T_ox, T_fu, V_air, V_fu, V_ox

  integer :: iN2  = -1
  integer :: iO2  = -1
  integer :: iH2  = -1
  integer :: iCH4 = -1
  integer :: iCO2 = -1
  integer :: iH2O = -1

contains

end module probdata_module
