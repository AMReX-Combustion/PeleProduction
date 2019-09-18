#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_3d_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none

  private

  public :: bcfunction, zero_visc, getZone, dcma_error

contains

  ! ::: -----------------------------------------------------------

  integer function getZone(x, y, z)

    use probdata_module, only : domnlo, BL_FUELPIPE, BL_COFLOW, BL_VOLUME

    REAL_T x, y, z

    getZone = BL_VOLUME
    if (z.gt.domnlo(3)) then
       getZone = BL_COFLOW
    else
       getZone = BL_FUELPIPE
    endif

  end function getZone

  ! ::: -----------------------------------------------------------

  subroutine bcfunction(x,y,z,dir,norm,time,u,v,w,rho,Yl,T,h,dx,getuv) &
       bind(C, name="bcfunction")

    use network,   only: nspecies
    use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_TfromHYpt
    use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, domnhi
    use probdata_module, only : blobr, bcinit, xfrontw, splitx, Tfrontw, &
         Y_bc, T_bc, u_bc, v_bc, w_bc
    use probdata_module, only : BL_FUELPIPE, BL_COFLOW, iN2, iO2, iNC12H26
    use fuego_chemistry, only : CKHMS

    REAL_T, intent(in) ::  x, y, z, time
    REAL_T, intent(inout) :: u, v, w, rho, Yl(0:*), T, h, dx(dim)
    logical, intent(in) :: getuv
    integer, intent(in) :: dir, norm

    REAL_T rho_temp(1), h_temp(1), T_temp(1)
    integer n, zone, len
    REAL_T eta, eta1, xmid, etamax, Patm
    REAL_T h_fu(0:nspecies-1), h_ox(0:nspecies-1), hmix
    REAL_T Wf, Wa, Wm, mf, Yf

    REAL_T,  parameter :: HtoTerrMAX = BL_REAL_E(7.8,-12)
    integer, parameter :: HtoTiterMAX = 20
    REAL_T res(0:HtoTiterMAX-1)
    integer Niter
    integer b(3)
    data  b / 1, 1, 1 /

    if (.not. bcinit) then
       call bl_abort('Need to initialize boundary condition function')
    end if

    eta = 0.5d0*(1.d0 - TANH(2.d0*(sqrt(x**2+y**2)-blobr)/Tfrontw))
    do n = 0, nspecies-1
       Yl(n) = Y_bc(n,BL_FUELPIPE)*eta + (1.d0-eta)*Y_bc(n,BL_COFLOW)
    end do
#if 0
    T = T_bc(BL_FUELPIPE)*eta + (1.d0-eta)*T_bc(BL_COFLOW)
#else

    call CKHMS(T_bc(BL_COFLOW), h_ox)
    h_fu = h_fu*1.d-4 ! cgs to MKS


    !   fuel enthalpy inncluding heat of vaporization
    h_fu( iNC12H26-1) = 10000.d0*((T_bc(BL_FUELPIPE)-298.d0)*0.375d0-352.1d0)/1.703348
    h_ox = h_ox*1.d-4 ! cgs to MKS

    hmix = Yl(iNC12H26-1)*h_fu(iNC12H26-1) &
         + Yl(iO2-1)     *h_ox(iO2-1) &
         + Yl(iN2-1)     *h_ox(iN2-1)

    T = T_bc(BL_COFLOW) ! initial guess
    call pphys_TfromHYpt(T,hmix,Yl,HtoTerrMAX,HtoTiterMAX,res,Niter)
#endif

    if (getuv .eqv. .TRUE.) then
       eta1 = 0.5d0*(1.d0 - TANH(2.d0*(sqrt(x**2+y**2)-splitx)/xfrontw))

       if(time.le.0.00006) then
          w = w_bc(BL_FUELPIPE)*eta1/2.0 + w_bc(BL_FUELPIPE)*eta1/2.0*time/0.0006
       else
          w = w_bc(BL_FUELPIPE)*eta1 + (1.d0-eta1)*w_bc(BL_COFLOW)
       endif
       u = u_bc(BL_FUELPIPE)*eta1 + (1.d0-eta1)*u_bc(BL_COFLOW)
       v = v_bc(BL_FUELPIPE)*eta1 + (1.d0-eta1)*v_bc(BL_COFLOW)
    endif

    Patm = pamb / pphys_getP1atm_MKS()
    T_temp(1) = T

    call pphys_RHOfromPTY(b, b, &
                          rho_temp(1), DIMARG(b), DIMARG(b), &
                          T_temp(1),   DIMARG(b), DIMARG(b), &
                          Yl, DIMARG(b), DIMARG(b), Patm)
    call pphys_HMIXfromTY(b, b, &
                          h_temp(1),   DIMARG(b), DIMARG(b), &
                          T_temp(1),   DIMARG(b), DIMARG(b), &
                          Yl, DIMARG(b), DIMARG(b))

    rho = rho_temp(1)
    h = h_temp(1)
    T = T_temp(1)

  end subroutine bcfunction

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

  subroutine dcma_error(tag, DIMS(tag), &
                        set, clear, hrr, DIMS(hrr), &
                        lo,hi,nvar, &
                        domlo, domhi, &
                        dx, xlo, &
                        problo,time, &
                        level, value) bind(C, name="dcma_error")
      use probdata_module
      implicit none
      integer   DIMDEC(tag)
      integer   DIMDEC(hrr)
      integer   nvar, set, clear, level
      integer   lo(DIM), hi(DIM)
      integer   domlo(DIM), domhi(DIM)
      REAL_T    dx(DIM), xlo(DIM), problo(DIM), time
      integer   tag(DIMV(tag))
      REAL_T    hrr(DIMV(hrr),nvar), value

      REAL_T x,y,z
      integer   i, j, k

      ! First component is mixfrac, second is HRR, third is OH, fourth is RO2
      do k = lo(3), hi(3)
        z = problo(3)+(k+0.5)*dx(3)
        do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            if(hrr(i,j,k,1).gt.mix_thresh.and.level.eq.0) tag(i,j,k)=set
            if(hrr(i,j,k,1).gt.mix_thresh.and.level.eq.1.and.abs(hrr(i,j,k,2)).gt.hrr_thresh) tag(i,j,k)=set

            if(level.eq.2) then
              if(z.lt.5*0.00017.and.hrr(i,j,k,1).gt.mix_thresh) then
                tag(i,j,k) = set
              endif

              if(z.lt.zmax_diff.and.abs(hrr(i,j,k,2)).gt.hrr_thresh) then
                tag(i,j,k) = set
              endif

              if(abs(hrr(i,j,k,2)).gt.hrr_thresh.and.hrr(i,j,k,3).lt.OH_thresh) tag(i,j,k)=set
              tag(i,j,k) = merge(set,tag(i,j,k), &
                                abs(hrr(i,j,k,4)).gt.RO2_thresh)

              if(z.gt.zmax_diff) then
                tag(i,j,k) = merge(clear,tag(i,j,k), &
                                  abs(hrr(i,j,k,3)).ge.OH_thresh)
              endif
            endif

            if(level.eq.3.and.hrr(i,j,k,4).ge.RO2_thresh) tag(i,j,k)=set
            if(level.eq.3.and.hrr(i,j,k,3).ge.OH_thresh.and.z.lt.0.002) tag(i,j,k)=set

            if(z.gt.zmax_mix) tag(i,j,k)=clear

          enddo
        enddo
      enddo
  end subroutine dcma_error

end module user_defined_fcts_3d_module
