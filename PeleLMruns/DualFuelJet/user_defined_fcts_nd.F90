#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include "mechanism.h"

module user_defined_fcts_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private

  public :: bcfunction, zero_visc!, dcma_error

  integer, parameter :: nspecies = NUM_SPECIES
  
contains

  ! ::: -----------------------------------------------------------

  subroutine bcfunction(x, dx, dir, norm, time, getuvw, &
                        vel, rho, Yl, T,  h) &
                        bind(C, name="bcfunction")

    use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor
    use mod_Fvar_def, only : pamb
    use probdata_module, only : BL_FUELPIPE, BL_COFLOW, &
                                bcinit, blobr, Tfrontw, splitx, xfrontw, &
                                Y_bc, T_bc, h_bc, rho_bc, u_bc, v_bc, w_bc, &
                                get_T_from_hY

    use PeleLM_nD, only: pphys_RHOfromPTY, pphys_TfromHY
    use PeleLM_F,  only: pphys_getP1atm_MKS
    
    REAL_T, intent(in)  :: x(3)
    REAL_T, intent(in)  :: dx(3)
    integer, intent(in) :: dir, norm  ! This specify the direction and orientation of the face
    REAL_T, intent(in)  :: time
    logical, intent(in) :: getuvw
    REAL_T, intent(out) :: vel(3)
    REAL_T, intent(out) :: rho
    REAL_T, intent(out) :: Yl(0:*)
    REAL_T, intent(out) :: T
    REAL_T, intent(out) :: h

    integer :: i
    REAL_T  :: eta, eta1
    logical :: is_ioproc

    REAL_T Patm

    integer :: Niter
    REAL_T,parameter :: Ttol = 1.d-12
    integer, parameter :: NiterMAX = 20
    REAL_T :: Tres(NiterMAX)
    integer,parameter,dimension(3) :: b_lo = (/1,1,1/)
    integer,parameter,dimension(3) :: b_hi = (/1,1,1/)
    REAL_T rarr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1))
    REAL_T harr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1))
    REAL_T Tarr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1))
    REAL_T Yarr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1),nspecies)


    Patm = pamb / pphys_getP1atm_MKS()
    
    if (.not. bcinit) then
       call bl_abort('Need to initialize boundary condition function')
    end if

    eta = 0.5d0*(1.d0 - TANH(2.d0*(sqrt(x(1)*x(1)+x(2)*x(2))-blobr)/Tfrontw))
    do i = 0, nspecies-1
       Yl(i) = Y_bc(i,BL_FUELPIPE)*eta + (1.d0-eta)*Y_bc(i,BL_COFLOW)
    end do
    h = h_bc(BL_FUELPIPE)*eta + (1.d0-eta)*h_bc(BL_COFLOW)
    T = T_bc(BL_FUELPIPE)*eta + (1.d0-eta)*T_bc(BL_COFLOW)

    Tarr(1,1,1) = T
    Yarr(1,1,1,1:nspecies) = Yl(0:nspecies-1)
    harr(1,1,1) = h
    Niter = pphys_TfromHY(b_lo, b_hi, &
                          Tarr, b_lo, b_hi, &
                          harr, b_lo, b_hi, &
                          Yarr, b_lo, b_hi, &
                          Ttol, NiterMAX, Tres)
    T = Tarr(1,1,1)

    if (getuvw) then
       eta1 = 0.5d0*(1.d0 - TANH(2.d0*(sqrt(x(1)*x(1)+x(2)*x(2))-splitx)/xfrontw))
       vel(1) = u_bc(BL_FUELPIPE)*eta1 + (1.d0-eta1)*u_bc(BL_COFLOW)
       vel(2) = v_bc(BL_FUELPIPE)*eta1 + (1.d0-eta1)*v_bc(BL_COFLOW)
       vel(3) = w_bc(BL_FUELPIPE)*eta1 + (1.d0-eta1)*w_bc(BL_COFLOW)
    endif

    
    Tarr(1,1,1) = T
    Yarr(1,1,1,1:nspecies) = Yl(0:nspecies-1)
    call pphys_RHOfromPTY(b_lo, b_hi, &
                          rarr, b_lo, b_hi, &
                          Tarr, b_lo, b_hi, &
                          Yarr, b_lo, b_hi, Patm)
    rho = rarr(1,1,1)

    ! --------------------------------------------------------------------------
    ! DEBUG Output
    ! --------------------------------------------------------------------------
!     is_ioproc = amrex_pd_ioprocessor()
!     if (is_ioproc) then
!       if (T .lt. 0.95*T_bc(BL_COFLOW)) then
!         write(6,'(2x,a)') 'Inlet boundary conditions:'
!         write(6,'(4x,a22,1x,f12.7)') 'u = ', u
!         write(6,'(4x,a22,1x,f12.7)') 'v = ', v
!         write(6,'(4x,a22,1x,f12.7)') 'w = ', w
!         write(6,'(4x,a22,1x,f12.7)') 'T = ', T
!         write(6,'(4x,a22,1x,f12.1)') 'h = ', h
!         write(6,'(4x,a22,1x,f12.7)') 'rho = ', rho
!         do i = 1, nspecies
!           if (Yl(i-1) .gt. 1e-14) then
!             write(6,'(4x,a22,1x,f12.7)') adjustl(spec_names(i)), Yl(i-1)
!           endif
!         enddo
!       endif
!     endif

  end subroutine bcfunction

  ! ::: -----------------------------------------------------------

  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
       dx,problo,bc,idir,isrz,id,ncomp) &
       bind(C, name="zero_visc")

    use mod_Fvar_def, only : Temp, FirstSpec, RhoH, LastSpec
    use probdata_module, only : domnlo

    implicit none
    integer DIMDEC(diff)
    integer lo(dim), hi(dim)
    integer domlo(dim), domhi(dim)
    integer bc(2*dim)
    integer idir, isrz, id, ncomp
    REAL_T  diff(DIMV(diff),*)
    REAL_T  dx(dim)
    REAL_T  problo(dim)

  end subroutine zero_visc

  ! ::: -----------------------------------------------------------


end module user_defined_fcts_nd_module
