#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>

module prob_nd_module

  use fuego_chemistry
  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort

  implicit none

  private
  
  public :: amrex_probinit, setupbc, init_data

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
    use mod_Fvar_def, only : pamb, fuelID, domnhi, domnlo
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

!------------------------------------

  subroutine setupbc() bind(c, name='setupbc')
      
    use network,   only: nspecies, spec_names
    use PeleLM_F, only: pphys_getP1atm_MKS
    use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo
    use probdata_module, only : Y_bc, T_bc, u_bc, v_bc, w_bc, rho_bc, h_bc
    use probdata_module, only : bcinit, iN2, iO2, iCO2, iCH4, iH2O, &
         BL_FUELPIPE, BL_OUTFLOW, BL_OXIDIZER, BL_AIR, BL_PIPEEND, BL_VOLUME,&
         T_air, T_ox, T_fu, V_air, V_fu, V_ox

    implicit none
    
    REAL_T Patm, pmf_vals(nspecies+3), a
    REAL_T Xt(nspecies), Yt(nspecies), loc
    REAL_T rho  ! added
    integer zone, n, fuelZone, airZone, oxZone, volZone, peZone, ofZone, region, len
    integer num_zones_defined, i
    integer :: b_lo(3), b_hi(3)
    data  b_lo(:) / 1, 1, 1 /
    data  b_hi(:) / 1, 1, 1 /

    Patm = pamb / pphys_getP1atm_MKS()

    do i = 1,nspecies
       if (spec_names(i) .eq. "O2") iO2 = i
       if (spec_names(i) .eq. "CH4") iCH4 = i
       if (spec_names(i) .eq. "N2") iN2 = i
       if (spec_names(i) .eq. "H2O") iH2O = i
       if (spec_names(i) .eq. "CO2") iCO2 = i
    enddo
    if (iO2.le.0 .or. iO2.gt.nspecies .or. &
         iCH4.le.0 .or. iCH4.gt.nspecies .or. &
         iN2.le.0 .or. iN2.gt.nspecies .or. &
         iH2O.le.0 .or. iH2O.gt.nspecies .or. &
         iCO2.le.0 .or. iCO2.gt.nspecies) then
       call bl_pd_abort('Could not find all special species in mech')
    endif

    !     A diffusion flame

    fuelZone = BL_FUELPIPE
    oxZone   = BL_OXIDIZER
    airZone  = BL_AIR
    volZone  = BL_VOLUME
    peZone   = BL_PIPEEND
    ofZone   = BL_OUTFLOW
    num_zones_defined = 6

    !     Fuel
    do n = 1,nspecies
       Yt(n) = 0.d0
    end do

    Yt(iCH4) = 1.0d0
    Yt(iO2) = 0.d0
    Yt(iH2O) = 0.d0
    Yt(iN2) = 1.d0 -  Yt(iCH4) -  Yt(iO2) - Yt(iH2O)

    do n=1,nspecies
       Y_bc(n-1,fuelZone) = Yt(n)
    end do
    T_bc(fuelZone) = T_fu
    u_bc(fuelZone) = 0.d0
    v_bc(fuelZone) = V_fu
    w_bc(fuelZone) = 0.d0

    !     Oxidizer
    do n = 1,nspecies
       Yt(n) = 0.d0
    end do

    Yt(iO2) = 0.2395d0
    Yt(iCO2) = 0.000d0
    Yt(iN2) = 1.d0 -  Yt(iO2) -  Yt(iCO2)
       
    do n=1,nspecies
       Y_bc(n-1,oxZone) = Yt(n)
    end do
    T_bc(oxZone) = T_ox
    u_bc(oxZone) = 0.d0
    v_bc(oxZone) = V_ox
    w_bc(oxZone) = 0.d0

    !     Air
    do n=1,nspecies
       Xt(n) = zero
    enddo
    Xt(iO2) = 0.21d0
    Xt(iN2) = 1.d0 - Xt(iO2)

    CALL CKXTY (Xt, Yt)         
    do n=1,nspecies
       Y_bc(n-1,airZone) = Yt(n)
    end do
       
    T_bc(airZone) = T_air
    u_bc(airZone) = 0.d0
    v_bc(airZone) = V_air
    w_bc(airZone) = 0.d0

    !     Pipeend (as oxidizer
    do n=1,nspecies
       Y_bc(n-1,peZone) = Y_bc(n-1,oxZone)
    end do
         
    T_bc(peZone) = T_bc(oxZone)
    u_bc(peZone) = u_bc(oxZone)
    v_bc(peZone) = 0.d0
    w_bc(peZone) = 0.d0

    !     Volume (as air)
    do n=1,nspecies
       Y_bc(n-1,volZone) = Y_bc(n-1,airZone)
    end do

    T_bc(volZone) = T_bc(airZone)
    u_bc(volZone) = u_bc(airZone)
    v_bc(volZone) = v_bc(airZone)
    w_bc(volZone) = 0.d0

    !     Outflow (as air, except phi)
    do n=1,nspecies
       Y_bc(n-1,ofZone) = Y_bc(n-1,airZone)
    end do
    
    T_bc(ofZone) = T_bc(airZone)
    u_bc(ofZone) = u_bc(airZone)
    v_bc(ofZone) = v_bc(airZone)
    w_bc(ofZone) = 0.d0

    do zone=1,num_zones_defined

!     Set density and hmix consistent with data

       call pphys_RHOfromPTY(b_lo, b_hi, &
                             rho_bc(zone), b_lo, b_hi,&
                             T_bc(zone),   b_lo, b_hi,&
                             Y_bc(0,zone), b_lo, b_hi, Patm)

       call pphys_HMIXfromTY(b_lo, b_hi, &
                             h_bc(zone),   b_lo, b_hi,&
                             T_bc(zone),   b_lo, b_hi,&
                             Y_bc(0,zone), b_lo, b_hi)

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
  subroutine init_data(level, time, lo, hi, nscal,&
                       vel, scal, s_lo, s_hi, &
                       press, p_lo, p_hi,&
                       delta, xlo, xhi) bind(c, name='init_data')

    use network,   only: nspecies
    use PeleLM_F,  only: pphys_getP1atm_MKS
    use PeleLM_nD, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac
    use mod_Fvar_def, only : bathID, domnhi, domnlo
    use probdata_module, only : BL_XLO, blobr, blobT, blobw, blobx, bloby
    use probdata_module, only : bcinit, iN2, iO2, iCO2, iCH4, iH2O, &
         BL_FUELPIPE, BL_OUTFLOW, BL_OXIDIZER, BL_AIR, BL_PIPEEND, BL_VOLUME
    use user_defined_fcts_nd_module, only : bcfunction

    implicit none

! In/Out
    integer, intent(in) :: level, nscal
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: p_lo(3), p_hi(3)
    REAL_T, intent(in)  :: xlo(3), xhi(3)
    REAL_T, intent(in)  :: time, delta(3)
    REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),dim), intent(out) :: vel
    REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal), intent(out) :: scal
    REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)), intent(out) :: press

! Local
    integer :: tmpi, nPMF
    integer :: i, j, k, n, airZone, fuelZone, zone, len
    REAL_T :: x, y, z, r, Yl(nspecies), Xl(nspecies), Patm
    REAL_T :: pmf_vals(nspecies+3), y1, y2, dx
    REAL_T :: pert,Lx,FORT_P1ATMMKS,eta,u,v,rho,T,h
    REAL_T, dimension(3) :: coor, velocity

    if (iN2.lt.1 .or. iN2.gt.nspecies) then
       call bl_pd_abort('N2 id not sest prior to calling INITDATA')
    endif
    
    do k = lo(3), hi(3)
       z = (float(k)+.5d0)*delta(3)+domnlo(3)
       do j = lo(2), hi(2)
          y = (float(j)+.5)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5)*delta(1)+domnlo(1)

             coor(1) = x  
             coor(2) = y  
             coor(3) = z  
             call bcfunction(coor, delta, 1,-1,time,.true.,velocity,rho,Yl,T,h)
             do n=1,nspecies
                scal(i,j,k,FirstSpec-1+n) = Yl(n)
             enddo
             scal(i,j,k,Temp) = T

             if (blobr.gt.0.d0) then
                r = SQRT((x-blobx)**2 + (y-bloby)**2)
                eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/blobw))
                scal(i,j,k,Temp) = blobT*eta + scal(i,j,k,Temp)*(1.d0-eta)
             endif

             vel(i,j,k,1) = velocity(1)
             vel(i,j,k,2) = velocity(2)
#if ( AMREX_SPACEDIM == 3 ) 
             vel(i,j,k,3) = velocity(3)
#endif
          enddo
       enddo
    enddo

    Patm = pamb / pphys_getP1atm_MKS()

    call pphys_RHOfromPTY(lo,hi, &
                          scal(:,:,:,Density),   s_lo, s_hi, &
                          scal(:,:,:,Temp),      s_lo, s_hi, &
                          scal(:,:,:,FirstSpec), s_lo, s_hi, &
                          Patm)
    call pphys_HMIXfromTY(lo,hi, &
                          scal(:,:,:,RhoH),      s_lo, s_hi, &
                          scal(:,:,:,Temp),      s_lo, s_hi, &
                          scal(:,:,:,FirstSpec), s_lo, s_hi)


    do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,nspecies-1
               scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n)*scal(i,j,k,Density)
            enddo
            scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
         enddo
      enddo
    enddo

  end subroutine init_data

end module prob_nd_module
