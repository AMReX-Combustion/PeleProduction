#include <AMReX_REAL.H>


module probdata_module

  use network, only : nspecies

  implicit none

  REAL_T :: standoff, pert_scale, pertmag, blobz
  REAL_T domnlo(3), domnhi(3)

  logical :: bcinit = .false.
  REAL_T :: u_bc, v_bc, w_bc
  REAL_T :: Y_bc(0:nspecies-1), T_bc(1), h_bc(1), rho_bc(1)
    
contains

end module probdata_module
