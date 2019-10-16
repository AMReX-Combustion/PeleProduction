#include <AMReX_REAL.H>


module probdata_module

  use network, only : nspecies

  implicit none

  ! from probdata.H
    REAL_T :: standoff
    REAL_T :: pertmag
    
    ! from bc.H
    
    logical :: bcinit
    
    REAL_T :: u_bc, v_bc, w_bc
    REAL_T :: Y_bc(0:nspecies-1,2), T_bc(1), h_bc(1,2), rho_bc(1,2)

    REAL_T :: midtanh
    REAL_T :: splitx
    REAL_T :: widthtanh
    REAL_T :: H2_enrich       ! H2 enrichment in volume (or mole fraction)
    REAL_T :: T_in
    REAL_T :: Zst
    
    integer, parameter :: flame_dir = 2
  
contains

!subroutines here

end module probdata_module
