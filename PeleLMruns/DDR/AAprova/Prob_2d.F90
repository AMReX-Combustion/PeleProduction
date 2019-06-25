#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>

module prob_2D_module

  use fuego_chemistry

  implicit none

  private
  
  public :: amrex_probinit, setupbc, init_data, getZone

contains

! ::: -----------------------------------------------------------
! ::: This routine is called at problem initialization time
! ::: and when restarting from a checkpoint file.
! ::: The purpose is (1) to specify the initial time value
! ::: (not all problems start at time=0.0) and (2) to read
! ::: problem specific data from a namelist or other input
! ::: files and possibly store them or derived information
! ::: in FORTRAN common blocks for later use.
! ::: 
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: init      => TRUE if called at start of problem run
! :::              FALSE if called from restart
! ::: strttime <=  start problem with this time variable
! ::: 
! ::: -----------------------------------------------------------
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c, name='amrex_probinit')

    use PeleLM_F,  only: pphys_getP1atm_MKS
    use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
    use mod_Fvar_def, only : fuelID, domnhi, domnlo, dim
    use probdata_module, only : rho_bc, Y_bc
    use probdata_module, only : T_fu, T_ox, T_air, pipeTh, pipeBL,&
         V_fu, V_ox, V_air, fuel_ox_split, ox_air_split,blobx, bloby, blobr, blobT, blobw

    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(dim), probhi(dim)

    integer i, istemp
    REAL_T FORT_P1ATMMKS, area

    namelist /fortin/ T_fu, T_ox, T_air, pipeTh, pipeBL,&
                      V_fu, V_ox, V_air, fuel_ox_split, ox_air_split,&
                      blobx, bloby, blobr, blobT, blobw
    namelist /heattransin/ pamb

!
!      Build 'probin' filename -- the name of file containing fortin namelist.
!
    integer maxlen, isioproc
    parameter (maxlen=256)
    character probin*(maxlen)

    call bl_pd_is_ioproc(isioproc)

    if (namlen .gt. maxlen) then
       call bl_pd_abort('probin file name too long')
    end if

    if (namlen .eq. 0) then
       namlen = 6
       probin(1:namlen) = 'probin'
    else
       do i = 1, namlen
          probin(i:i) = char(name(i))
       end do
    endif

!     Load domain dimensions into common (put something computable there for SDIM<3)
    domnlo(:) = 0.d0
    domnhi(:) = 0.d0

    domnlo(1:dim) = problo(1:dim)
    domnhi(1:dim) = probhi(1:dim)

    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')
      
    ! Set defaults
    pamb = pphys_getP1atm_MKS()

    blobr = 0.d0
    blobx = 0.5d0 * (domnhi(1) + domnlo(1))
    bloby = 0.5d0 * (domnhi(2) + domnlo(2))
    blobw = 1.d0
    fuel_ox_split = .001d0
    ox_air_split =  .0125d0
    pipeTh = 0.d0
    pipeBL = 1.d0
    V_fu = -0.2d0
    V_ox = 0.2d0
    V_air = 0.2d0
    T_fu = 300.d0
    T_ox = 300.d0
    T_air = 300.d0

    ! Note: for setup with no coflow, set Ro=Rf+wallth

    if (isioproc .eq. 1) write(6,*)"reading fortin"
    read(untin,fortin)
    if (isioproc .eq. 1) write(6,*)"done reading fortin"

    read(untin,heattransin)
 
!     Set up boundary functions
    if (isioproc .eq. 1) write(6,*)" setup bc"
    call setupbc()
    if (isioproc .eq. 1) write(6,*)" done setup bc"

    area = 1.d0
    do i=1,dim-1
       area = area*(domnhi(i)-domnlo(i))
    enddo

    if (isioproc.eq.1) then
       write(6,fortin)
       write(6,heattransin)
#if defined(BL_DO_FLCT)
       write(6,flctin)
#endif
    end if

  end subroutine amrex_probinit

  subroutine setupbc() bind(c, name='setupbc')
      
    use network,   only: nspec
    use PeleLM_F, only: pphys_getP1atm_MKS, pphys_get_spec_name2
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, maxspec, maxspnml, dim
    use probdata_module, only : Y_bc, T_bc, u_bc, v_bc, rho_bc, h_bc
    use probdata_module, only : bcinit, iN2, iO2, iCO2, iCH4, iH2O, &
         BL_FUELPIPE, BL_OUTFLOW, BL_OXIDIZER, BL_AIR, BL_PIPEEND, BL_VOLUME,&
         T_air, T_ox, T_fu, V_air, V_fu, V_ox
    implicit none
    
    REAL_T Patm, pmf_vals(maxspec+3), a
    REAL_T Xt(maxspec), Yt(maxspec), loc
    REAL_T rho	! added
    integer zone, n, getZone, fuelZone, airZone, oxZone, volZone, peZone, ofZone, region, len
    integer b(dim)
    integer num_zones_defined
    character*(maxspnml) name

    b = 1
    Patm = pamb / pphys_getP1atm_MKS()

    do n=1,Nspec
       call pphys_get_spec_name2(name,n)
       if (name .eq. 'N2' ) iN2 = n
       if (name .eq. 'O2' ) iO2 = n
       if (name .eq. 'CO2' ) iCO2 = n
       if (name .eq. 'CH4' ) iCH4 = n
       if (name .eq. 'H2O' ) iH2O = n
    enddo

    !     A diffusion flame

    fuelZone = BL_FUELPIPE
    oxZone   = BL_OXIDIZER
    airZone  = BL_AIR
    volZone  = BL_VOLUME
    peZone   = BL_PIPEEND
    ofZone   = BL_OUTFLOW
    num_zones_defined = 6

    !     Fuel
    do n = 1,Nspec
       Yt(n) = 0.d0
    end do

    Yt(iCH4) = 1.0d0
    Yt(iO2) = 0.d0
    Yt(iH2O) = 0.d0
    Yt(iN2) = 1.d0 -  Yt(iCH4) -  Yt(iO2) - Yt(iH2O)

    do n=1,Nspec
       Y_bc(n-1,fuelZone) = Yt(n)
    end do
    T_bc(fuelZone) = T_fu
    u_bc(fuelZone) = 0.d0
    v_bc(fuelZone) = V_fu

    !     Oxidizer
    do n = 1,Nspec
       Yt(n) = 0.d0
    end do

    Yt(iO2) = 0.2395d0
    Yt(iCO2) = 0.000d0
    Yt(iN2) = 1.d0 -  Yt(iO2) -  Yt(iCO2)
       
    do n=1,Nspec
       Y_bc(n-1,oxZone) = Yt(n)
    end do
    T_bc(oxZone) = T_ox
    u_bc(oxZone) = 0.d0
    v_bc(oxZone) = V_ox

    !     Air
    do n=1,Nspec
       Xt(n) = zero
    enddo
    Xt(iO2) = 0.21d0
    Xt(iN2) = 1.d0 - Xt(iO2)

    CALL CKXTY (Xt, Yt)         
    do n=1,Nspec
       Y_bc(n-1,airZone) = Yt(n)
    end do
       
    T_bc(airZone) = T_air
    u_bc(airZone) = 0.d0
    v_bc(airZone) = V_air

    !     Pipeend (as oxidizer
    do n=1,Nspec
       Y_bc(n-1,peZone) = Y_bc(n-1,oxZone)
    end do
         
    T_bc(peZone) = T_bc(oxZone)
    u_bc(peZone) = u_bc(oxZone)
    v_bc(peZone) = 0.d0

    !     Volume (as air)
    do n=1,Nspec
       Y_bc(n-1,volZone) = Y_bc(n-1,airZone)
    end do

    T_bc(volZone) = T_bc(airZone)
    u_bc(volZone) = u_bc(airZone)
    v_bc(volZone) = v_bc(airZone)

    !     Outflow (as air, except phi)
    do n=1,Nspec
       Y_bc(n-1,ofZone) = Y_bc(n-1,airZone)
    end do
    
    T_bc(ofZone) = T_bc(airZone)
    u_bc(ofZone) = u_bc(airZone)
    v_bc(ofZone) = v_bc(airZone)

    do zone=1,num_zones_defined

!     Set density and hmix consistent with data

       call pphys_RHOfromPTY(b, b, &
                             rho_bc(zone), DIMARG(b), DIMARG(b),&
                             T_bc(zone),   DIMARG(b), DIMARG(b),&
                             Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)

       call pphys_HMIXfromTY(b, b, &
                             h_bc(zone),   DIMARG(b), DIMARG(b),&
                             T_bc(zone),   DIMARG(b), DIMARG(b),&
                             Y_bc(0,zone), DIMARG(b), DIMARG(b))

    enddo

    bcinit = .true.

  end subroutine setupbc


! ::: -----------------------------------------------------------
! ::: This routine is called at problem setup time and is used
! ::: to initialize data on each grid.  The velocity field you
! ::: provide does not have to be divergence free and the pressure
! ::: field need not be set.  A subsequent projection iteration
! ::: will define aa divergence free velocity field along with a
! ::: consistant pressure.
! ::: 
! ::: NOTE:  all arrays have one cell of ghost zones surrounding
! :::        the grid interior.  Values in these cells need not
! :::        be set here.
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: level     => amr level of grid
! ::: time      => time at which to init data             
! ::: lo,hi     => index limits of grid interior (cell centered)
! ::: nscal     => number of scalar quantities.  You should know
! :::		   this already!
! ::: vel      <=  Velocity array
! ::: scal     <=  Scalar array
! ::: press    <=  Pressure array
! ::: delta     => cell size
! ::: xlo,xhi   => physical locations of lower left and upper
! :::              right hand corner of grid.  (does not include
! :::		   ghost region).
! ::: -----------------------------------------------------------
  subroutine init_data(level,time,lo,hi,nscal,&
                       vel,scal,DIMS(state),press,DIMS(press),&
                       delta,xlo,xhi) bind(c, name='init_data')
    use network,   only: nspec
    use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
    use mod_Fvar_def, only : bathID, domnhi, domnlo, maxspec, maxspnml
    use probdata_module, only : BL_XLO, blobr, blobT, blobw, blobx, bloby
    use probdata_module, only : bcinit, iN2, iO2, iCO2, iCH4, iH2O, &
         BL_FUELPIPE, BL_OUTFLOW, BL_OXIDIZER, BL_AIR, BL_PIPEEND, BL_VOLUME
    use user_defined_fcts_2d_module, only : bcfunction

    implicit none

    integer    level, nscal
    integer    lo(dim), hi(dim)
    integer    DIMDEC(state)
    integer    DIMDEC(press)
    REAL_T     xlo(dim), xhi(dim)
    REAL_T     time, delta(dim)
    REAL_T     vel(DIMV(state),dim)
    REAL_T    scal(DIMV(state),nscal)
    REAL_T   press(DIMV(press))
    integer tmpi, nPMF

    integer i, j, n, airZone, fuelZone, getZone, zone, len
    REAL_T x, y, r, Yl(maxspec), Xl(maxspec), Patm
    REAL_T pmf_vals(maxspec+3), y1, y2, dx
    REAL_T pert,Lx,FORT_P1ATMMKS,eta,u,v,rho,T,h

    if (iN2.lt.1 .or. iN2.gt.Nspec) then
       call bl_pd_abort('N2 id not sest prior to calling INITDATA')
    endif
    
    do j = lo(2), hi(2)
       y = (float(j)+.5)*delta(2)+domnlo(2)
       do i = lo(1), hi(1)
          x = (float(i)+.5)*delta(1)+domnlo(1)

          call bcfunction(x,y,1,-1,time,u,v,rho,Yl,T,h,delta,.true.)

          do n=1,Nspec
             scal(i,j,FirstSpec-1+n) = Yl(n)
          enddo
          scal(i,j,Temp) = T

          if (blobr.gt.0.d0) then
             r = SQRT((x-blobx)**2 + (y-bloby)**2)
             eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/blobw))
             scal(i,j,Temp) = blobT*eta + scal(i,j,Temp)*(1.d0-eta)
          endif

          vel(i,j,1) = u
          vel(i,j,2) = v
          scal(i,j,Trac) = 0.d0

       enddo
    enddo

    Patm = pamb / pphys_getP1atm_MKS()

    call pphys_RHOfromPTY(lo,hi,&
                          scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state),&
                          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),&
                          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state),&
                          Patm)

    call pphys_HMIXfromTY(lo,hi,&
                          scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state),&
                          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state),&
                          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state))

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          do n = 0,Nspec-1
             scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
          enddo
          scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
       enddo
    enddo
  end subroutine init_data

  integer function getZone(x, y)
    
    use mod_Fvar_def, only : domnhi, domnlo, dim
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
      
end module prob_2D_module
