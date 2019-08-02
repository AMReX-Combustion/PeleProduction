#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_2d_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none
  
  private
  
  public :: getZone, bcfunction, zero_visc

contains

  integer function getZone(x, y)
    
    use mod_Fvar_def, only : domnhi, domnlo
    use probdata_module, only : BL_FUELPIPE, BL_OUTFLOW, BL_OXIDIZER, BL_AIR, BL_PIPEEND, BL_VOLUME,&
         fuel_ox_split, ox_air_split, pipeTh

    implicit none
    REAL_T x, y
    getZone = BL_VOLUME
    if (y.le.domnlo(2)) then
       if (ABS(x).le.fuel_ox_split) then
          getZone = BL_FUELPIPE
       else if (ABS(x).le.fuel_ox_split+pipeTh) then
          getZone = BL_PIPEEND
       else if (ABS(x).le.ox_air_split) then
          getZone = BL_OXIDIZER
       else
          getZone = BL_AIR
       endif
    else if (y.ge.domnhi(2)) then
       getZone = BL_OUTFLOW
    endif
  end function getZone
      
  subroutine bcfunction(x,y,coord,lohi,time,u,v,rho,Yl,T,h,dx,getuv) bind(c, name='bcfunction')

    use network,   only: nspecies
    use PeleLM_F,  only: pphys_getP1atm_MKS
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb
    use probdata_module, only : bcinit, rho_bc, Y_bc, T_bc, h_bc, v_bc, &
         BL_FUELPIPE, BL_OUTFLOW, BL_OXIDIZER, BL_AIR, BL_PIPEEND, BL_VOLUME,&
         fuel_ox_split, blobw, bloby, pipeBL, pipeTh

    integer coord, lohi
    REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(2), r
    logical getuv
    integer b(2)

    integer n, zone
    REAL_T eta, xmid, Patm, rhoV(1,1), TV(1,1), hV(1,1)
         
    REAL_T Wf, Wa, Wm, mf, Yf

    REAL_T,  parameter :: HtoTerrMAX = BL_REAL_E(7.8,-12)
    integer, parameter :: HtoTiterMAX = 20
    REAL_T res(0:HtoTiterMAX-1)
    integer Niter
    integer zoneL, zoneR, zoneT

    data  b / 1, 1 /
    
    if (.not. bcinit) then
       call bl_abort('Need to initialize boundary condition function')
    end if

    zone = getZone(x,y)
    if (zone .eq. BL_OUTFLOW) then
       rho = rho_bc(zone)
       do n = 0, nspecies-1
          Yl(n) = Y_bc(n,zone)
       end do
       T = T_bc(zone)
       h = h_bc(zone)
       if (getuv .eqv. .TRUE.) then
          u = zero
          v = v_bc(zone)
       endif
    else
       zoneL = BL_FUELPIPE
       zoneR = BL_OXIDIZER
       zoneT = BL_VOLUME         
       eta = 0.5d0*(1.d0 - TANH(2.d0*(ABS(x)-fuel_ox_split)/blobw))
       do n = 0, nspecies-1
          Yl(n) = Y_bc(n,zoneL)*eta + Y_bc(n,zoneR)*(1.d0-eta)
       end do
       T = T_bc(zoneL)*eta + T_bc(zoneR)*(1.d0-eta)

       eta = 0.5d0*(1.d0 - TANH(2.d0*(y-bloby)/blobw))
       do n = 0, nspecies-1
          Yl(n) = Yl(n)*eta + Y_bc(n,zoneT)*(1.d0-eta)
       end do
       T = T*eta + T_bc(zoneT)*(1.d0-eta)

       Patm = pamb / pphys_getP1atm_MKS()

       TV(b(1),b(2)) = T
       call pphys_RHOfromPTY(b, b, &
                             rhoV, DIMARG(b), DIMARG(b),&
                             TV,   DIMARG(b), DIMARG(b),&
                             Yl,  DIMARG(b), DIMARG(b), Patm)
       rho = rhoV(b(1),b(2))
       
       call pphys_HMIXfromTY(b, b, &
                             hV,   DIMARG(b), DIMARG(b),&
                             TV,   DIMARG(b), DIMARG(b),&
                             Yl,  DIMARG(b), DIMARG(b))
       h = hV(b(1),b(2))
          
       if (getuv .eqv. .TRUE.) then
          u = zero
          zone = getZone(x,y)
          if (zone.eq.BL_FUELPIPE) then
             zoneL = BL_PIPEEND
             zoneR = BL_FUELPIPE
             r = fuel_ox_split - ABS(x)
	     eta = (1.d0 - TANH(2.d0*r/pipeBL))
	     v = v_bc(zoneR)*2.d0*(1-((ABS(x))**(2)/(fuel_ox_split**(2))))
!             eta = 0.5d0*(1.d0 - TANH(2.d0*r/pipeBL))
!             v = v_bc(zoneL)*eta + v_bc(zoneR)*(1.d0-eta)
             if (eta.lt.0.d0 .or. eta.gt.1.d0) then
                call bl_pd_abort()
             endif
          else if (zone.eq.BL_OXIDIZER) then
             zoneL = BL_PIPEEND
             zoneR = BL_OXIDIZER
             r = ABS(x) - (fuel_ox_split+pipeTh)
             eta = (1.d0 - TANH(2.d0*r/pipeBL))
!             eta = 0.5d0*(1.d0 - TANH(2.d0*r/pipeBL))
             v = v_bc(zoneL)*eta + v_bc(zoneR)*(1.d0-eta)
             if (eta.lt.0.d0 .or. eta.gt.1.d0) then
                call bl_pd_abort()
             endif
          else
             v = v_bc(zone)
          endif
       endif
    endif

  end subroutine bcfunction

! ::: -----------------------------------------------------------
! ::: This routine will zero out diffusivity on portions of the
! ::: boundary that are inflow, allowing that a "wall" block
! ::: the complement aperture
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: diff      <=> diffusivity on edges
! ::: DIMS(diff) => index extent of diff array
! ::: lo,hi      => region of interest, edge-based
! ::: domlo,hi   => index extent of problem domain, edge-based
! ::: dx         => cell spacing
! ::: problo     => phys loc of lower left corner of prob domain
! ::: bc         => boundary condition flag (on orient)
! :::                   in BC_TYPES::physicalBndryTypes
! ::: idir       => which face, 0=x, 1=y
! ::: isrz       => 1 if problem is r-z
! ::: id         => index of state, 0=u
! ::: ncomp      => components to modify
! :::
! ::: -----------------------------------------------------------

  subroutine zero_visc(diff,DIMS(diff),lo,hi,domlo,domhi, &
                           dx,problo,bc,idir,isrz,id,ncomp) &
                           bind(C, name="zero_visc")   
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, LastSpec
    use probdata_module, only : BL_PIPEEND
      
    implicit none
    integer DIMDEC(diff)
    integer lo(dim), hi(dim)
    integer domlo(dim), domhi(dim)
    integer bc(2*dim)
    integer idir, isrz, id, ncomp
    REAL_T  diff(DIMV(diff),*)
    REAL_T  dx(dim)
    REAL_T  problo(dim)

    integer i, j, n, Tid, RHid, YSid, YEid, ys, ye
    logical do_T, do_RH, do_Y
    REAL_T xl, xr, xh, y

    Tid  = Temp      - id + dim
    RHid = RhoH      - id + dim
    YSid = FirstSpec - id + dim
    YEid = LastSpec  - id + dim
         
    do_T  = (Tid  .GE. 1) .AND. (Tid  .LE. ncomp)
    do_RH = (RHid .GE. 1) .AND. (RHid .LE. ncomp)
    ys = MAX(YSid,1)
    ye = MIN(YEid,ncomp)
    do_Y = (ye - ys + 1) .GE. 1
    !     
    !     Do species, Temp, rhoH
    !     
    if ((idir.EQ.1) .AND. (lo(2) .LE. domlo(2))&
         .AND. (do_T .OR. do_RH .OR. do_Y) ) then
               
       y = float(j)*dx(2)+problo(2)
       j = lo(2)
       do i = lo(1), hi(1)
          
          xl = float(i)*dx(1)+problo(1) 
          xr = (float(i)+1.d0)*dx(1)+problo(1) 
          xh = 0.5d0*(xl+xr)
                  
          if ( (getZone(xl,y).eq.BL_PIPEEND) .OR.&
               (getZone(xh,y).eq.BL_PIPEEND) .OR.&
               (getZone(xr,y).eq.BL_PIPEEND)  ) then
                
             !if (do_T)  diff(i,j,Tid ) = 0.d0
             !if (do_RH) diff(i,j,RHid) = 0.d0
             if (do_Y) then
                do n=ys,ye
                   diff(i,j,n) = 0.d0
                enddo
             endif
                     
          endif
       end do
    endif

  end subroutine zero_visc
  
end module user_defined_fcts_2d_module

