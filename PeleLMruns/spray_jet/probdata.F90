#include <AMReX_REAL.H>

#include "mechanism.h"

module probdata_module

  implicit none

    ! from bc.H
    
    logical :: bcinit
    
    REAL_T :: u_bc, v_bc, w_bc
    REAL_T :: Y_bc(0:NUM_SPECIES-1), T_bc(1), h_bc(1), rho_bc(1)
    REAL_T :: X_O2_oxid, T_in, splitx, xfrontw
    REAL_T :: turb_scale = 1.d0
    
contains

!subroutines here

end module probdata_module
