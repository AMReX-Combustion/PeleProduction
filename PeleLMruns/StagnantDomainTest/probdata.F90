#include <AMReX_REAL.H>


module probdata_module

  use network, only : nspecies

  implicit none

    ! from probdata.H
    REAL_T :: T_in  
    
    ! from bc.H
    
    logical :: bcinit
    integer, parameter :: nzones = 1
    integer, parameter :: BL_AIR = 1
    
    REAL_T :: u_bc(nzones), v_bc(nzones), w_bc(nzones)
    REAL_T :: Y_bc(0:nspecies-1,nzones), T_bc(nzones), h_bc(nzones), rho_bc(nzones)
      
contains

!subroutines here

end module probdata_module
