
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>


module prob_2D_module

  use amrex_fort_module, only : dim=>amrex_spacedim
  use fuego_chemistry

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

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)
  
      use PeleLM_F,  only: pphys_getP1atm_MKS

      use mod_Fvar_def, only : pamb
      use mod_Fvar_def, only : fuelID, domnhi, domnlo

      use mod_Fvar_def, only : ac_hist_file, cfix, changemax_control, &
                               coft_old, controlvelmax, corr, dv_control, &
                               h_control, navg_pnts, scale_control, sest, &
                               tau_control, tbase_control, V_in, v_in_old, zbase_control, &
                               pseudo_gravity
      use probdata_module, only : rho_bc, Y_bc
      use probdata_module, only : T_in
      
      
      implicit none
      
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i,istemp
      REAL_T area

      namelist /fortin/ T_in

      namelist /heattransin/ pamb

      namelist /control/ tau_control, sest, cfix, changeMax_control, h_control, &
          zbase_control, pseudo_gravity, istemp,corr,controlVelMax,navg_pnts

!
!      Build `probin' filename -- the name of file containing fortin namelist.
!
      integer maxlen, isioproc
      parameter (maxlen=256)
      character probin*(maxlen)

      call bl_pd_is_ioproc(isioproc)

      if (init.ne.1) then
!         call bl_abort('probinit called with init ne 1')
      end if

      if (namlen .gt. maxlen) then
         call bl_abort('probin file name too long')
      end if

      if (namlen .eq. 0) then
         namlen = 6
         probin(1:namlen) = 'probin'
      else
         do i = 1, namlen
            probin(i:i) = char(name(i))
         end do
      endif

      untin = 9
      open(untin,file=probin(1:namlen),form='formatted',status='old')
      
!     Set defaults
      pamb = pphys_getP1atm_MKS()

      zbase_control = 0.d0

!     Initialize control variables
      tau_control = one
      sest = zero
      corr = one
      changeMax_control = .05
      coft_old = -one
      cfix = zero
      ac_hist_file = 'AC_History'
      h_control = -one
      dV_control = zero
      tbase_control = zero
      h_control = -one
      pseudo_gravity = 0
      istemp = 0
      navg_pnts = 10

      read(untin,fortin)
      
      read(untin,heattransin)
 
      read(untin,control)
      close(unit=untin)

!     Set up boundary functions
      call setupbc()
      

      if (isioproc.eq.1) then
         write(6,fortin)
         write(6,heattransin)
         write(6,control)
      end if

  end subroutine amrex_probinit
  
!------------------------------------
  
  subroutine setupbc()bind(C, name="setupbc")

    use network,   only: nspecies
    use PeleLM_F, only: pphys_getP1atm_MKS
    use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, V_in
    use mod_Fvar_def, only : fuelID, oxidID, bathID
    use probdata_module, only : Y_bc, T_bc, u_bc, v_bc, rho_bc, h_bc
    use probdata_module, only : bcinit
    use probdata_module, only : nzones, BL_AIR
    use probdata_module, only : T_in
  
    implicit none

    REAL_T Patm, pmf_vals(nspecies+3), Xt(nspecies), Yt(nspecies), loc
    
    integer n
    integer b(2)
    integer zone
    data  b / 1, 1 /
      
    Patm = pamb / pphys_getP1atm_MKS()

!   Ambient conditions
    do n = 1,nspecies
      Xt(n) = 0.d0
    enddo
    Xt(oxidID) = 0.21d0
    Xt(bathID) = 1.d0 - Xt(oxidID)

    CALL CKXTY (Xt, Yt)

    do n = 1,nspecies
      Y_bc(n-1,BL_AIR) = Yt(n)
    enddo
    T_bc(BL_AIR) = T_in
    u_bc(BL_AIR) = zero
    v_bc(BL_AIR) = zero
                           

!     Set density and hmix consistent with data
    do zone=1,nzones
      call pphys_RHOfromPTY(b, b, &
                           rho_bc(zone), DIMARG(b), DIMARG(b), &
                           T_bc(zone),   DIMARG(b), DIMARG(b), &
                           Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)
      call pphys_HMIXfromTY(b, b, &
                           h_bc(zone),   DIMARG(b), DIMARG(b), &
                           T_bc(zone),   DIMARG(b), DIMARG(b), &
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

  subroutine init_data(level,time,lo,hi,nscal, &
                       vel,scal,DIMS(state),press,DIMS(press), &
                       delta,xlo,xhi) &
                       bind(C, name="init_data")
                              

      use network,   only: nspecies
      use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb
      use mod_Fvar_def, only : bathID, domnhi, domnlo
      use probdata_module, only : BL_AIR, Y_bc
      use probdata_module, only : T_in

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
      integer nPMF

      integer i, j, n
      REAL_T x, y, Yl(nspecies), Xl(nspecies), Patm
      REAL_T pmf_vals(nspecies+3), y1, y2

!      write(6,*)" made it to initdata"
      if (bathID.lt.1 .or. bathID.gt.nspecies) then
         call bl_pd_abort()
      endif

      do j = lo(2), hi(2)
        y = (float(j)+.5d0)*delta(2)+domnlo(2)
        do i = lo(1), hi(1)
          x = (float(i)+.5d0)*delta(1)+domnlo(1)

          do n=1,nspecies
             scal(i,j,FirstSpec+n-1) = Y_bc(n-1,BL_AIR)
          enddo
          scal(i,j,Temp) = T_in

          vel(i,j,1) = zero
          vel(i,j,2) = zero

        end do
      end do

      Patm = pamb / pphys_getP1atm_MKS()
!      write(6,*)"Patm",Patm

      call pphys_RHOfromPTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),Density),  DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state), &
          Patm)

      call pphys_HMIXfromTY(lo,hi, &
          scal(ARG_L1(state),ARG_L2(state),RhoH),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),Temp),     DIMS(state), &
          scal(ARG_L1(state),ARG_L2(state),FirstSpec),DIMS(state)) 

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            do n = 0,nspecies-1
              scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n)*scal(i,j,Density)
            enddo
            scal(i,j,RhoH) = scal(i,j,RhoH)*scal(i,j,Density)
         enddo
      enddo
      
  end subroutine init_data
      


end module prob_2D_module
