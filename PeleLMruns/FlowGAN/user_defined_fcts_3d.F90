#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_3d_module

  implicit none
  
  private
  
  public :: bcfunction, zero_visc

contains
  
  subroutine bcfunction(x,y,z,dir,norm,time,u,v,w,rho,Yl,T,h,dx,getuvw) &
       bind(C, name="bcfunction")

    use mod_Fvar_def, only : dim
    use network,      only : nspec
    use probdata_module, only : u_bc, v_bc, w_bc, rho_bc, T_bc, h_bc, Y_bc
       
    implicit none

    REAL_T x, y, z, time, u, v, w, rho, Yl(0:*), T, h, dx(dim)
    integer dir, norm  ! This specify the direction and orientation of the face
    logical getuvw

    integer zone, n

    zone = 1
    rho = rho_bc(zone)
    do n = 0, nspec-1
       Yl(n) = Y_bc(n)
    end do
    T = T_bc(zone)
    h = h_bc(zone)

    if (getuvw .eqv. .TRUE.) then
       u = u_bc
       v = v_bc
       w = w_bc
    endif

  end subroutine bcfunction

  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
                           dx,problo,bc,idir,isrz,id,ncomp) &
                           bind(C, name="zero_visc")   

      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, LastSpec
      use mod_Fvar_def, only : domnhi, domnlo, dim
      
      implicit none
      integer DIMDEC(diff)
      integer lo(dim), hi(dim)
      integer domlo(dim), domhi(dim)
      integer bc(2*dim)
      integer idir, isrz, id, ncomp
      REAL_T  diff(DIMV(diff),*)
      REAL_T  dx(dim)
      REAL_T  problo(dim)



! Routine compiled but should be set by the user
! if there is a mix of inflox/wall at a boundary



  end subroutine zero_visc

end module user_defined_fcts_3d_module

