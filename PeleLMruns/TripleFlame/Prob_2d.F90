
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
      use probdata_module, only : standoff, pertmag, T_in, rho_bc, Y_bc, splitx
      use probdata_module, only : flame_dir, midtanh, widthtanh, H2_enrich
      use user_defined_fcts_2d_module, only: set_Zst
      
      
      implicit none
      
      integer init, namlen
      integer name(namlen)
      integer untin
      REAL_T problo(dim), probhi(dim)

      integer i,istemp
      REAL_T area

      namelist /fortin/ V_in, T_in, standoff, pertmag, &
                        midtanh, widthtanh, H2_enrich
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

!     Note: for setup with no coflow, set Ro=Rf+wallth
      standoff = zero
      pertmag = 0.d0
      H2_enrich = 0.0d0
      T_in = 300.0d0
      midtanh = 0.6*(domnhi(1)+domnlo(1))    ! Default is middomain in x
      widthtanh = 0.05*(domnhi(1)-domnlo(1))  ! Default is 1/10th of domain width in x 

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
      
!     Initialize control variables that depend on fortin variables
      V_in_old = V_in
      
      read(untin,heattransin)
 
      read(untin,control)
      close(unit=untin)

!     Set Zst
      call set_Zst()

!     Set up boundary functions
      call setupbc()

      splitx = 0.5d0 * (domnhi(1) + domnlo(1))
      
      area = 1.d0
      do i=1,dim
        if (flame_dir /= i) then
         area = area*(domnhi(i)-domnlo(i))
        endif
      enddo
      !scale_control = (Y_bc(fuelID-1,1)*rho_bc(1,1)+Y_bc(fuelID-1,2)*rho_bc(1,2)) * 0.5 * area
      scale_control = 1.0

      if (h_control .gt. zero) then
         !cfix = scale_control * h_control
         cfix = h_control
      endif

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
    use probdata_module, only : standoff, Y_bc, T_bc, u_bc, v_bc, rho_bc, h_bc
    use probdata_module, only : bcinit, T_in
    use user_defined_fcts_2d_module, only: set_Y_from_ksi
  
    implicit none

    REAL_T :: Patm
    
    integer n
    integer b(2)
    data  b / 1, 1 /
      
    Patm = pamb / pphys_getP1atm_MKS()
             
  ! Get the lean and rich mixture comp
    
    call set_Y_from_ksi(1.0d0,Y_bc(0:nspecies-1,1))     
    call set_Y_from_ksi(0.0d0,Y_bc(0:nspecies-1,2))     

    T_bc = T_in
    v_bc = V_in
    u_bc = zero
              
!   Set density and hmix consistent with data

    call pphys_RHOfromPTY(b, b, &
                         rho_bc(1,1), DIMARG(b), DIMARG(b), &
                         T_bc(1),     DIMARG(b), DIMARG(b), &
                         Y_bc(0,1),   DIMARG(b), DIMARG(b), Patm)
    call pphys_RHOfromPTY(b, b, &
                         rho_bc(1,2), DIMARG(b), DIMARG(b), &
                         T_bc(1),     DIMARG(b), DIMARG(b), &
                         Y_bc(0,2),   DIMARG(b), DIMARG(b), Patm)
    call pphys_HMIXfromTY(b, b, &
                         h_bc(1,1), DIMARG(b), DIMARG(b), &
                         T_bc(1),   DIMARG(b), DIMARG(b), &
                         Y_bc(0,1), DIMARG(b), DIMARG(b))
    call pphys_HMIXfromTY(b, b, &
                         h_bc(1,2), DIMARG(b), DIMARG(b), &
                         T_bc(1),   DIMARG(b), DIMARG(b), &
                         Y_bc(0,2), DIMARG(b), DIMARG(b))

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
      use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac
      use mod_Fvar_def, only : bathID, oxidID, domnhi, domnlo, V_in
      use probdata_module, only : Y_bc, rho_bc, T_bc, midtanh, widthtanh, splitx
      use user_defined_fcts_2d_module, only : set_Y_from_ksi

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
      REAL_T x, y, Yl(0:nspecies-1,2), Xl(nspecies), Patm
      REAL_T pmf_vals(nspecies+3), y1, y2
      REAL_T pert,Lx, tanhval, rad, y_lo
      REAL_T :: Gaussian_T, Gauss_T_width, Gauss_maxT
      REAL_T :: Gaussian_Spec, Gaussian_Spec_width
      REAL_T :: sumspec

!--------------------------------------------------------    
! Data initialization for triple flame case: 
!     - first set-up a mixing layer with a tanh in mixture fraction
!       between 0 and 1. It uses the fuel declared in input file. A
!       constant temperature is set first.
!     - Then a widening in the y-dir Gaussian profile of temperature
!       is super-imposed on the base flow. The composition in the hot
!       is replaced by air to avoid preheating premixed gas.
!         * y_lo : lower tip of the hot region
!         * Gauss_maxT : maximum temperature
!         * Gauss_T_width : width of the temperature profile
!         * Gaussian_Spec_width : width of the hot air profile
!
!     - Setting Gaussian_Spec_width > Gauss_T_width enable to avoid
!       heating fuel/air mixture at the border of the hot region.
!--------------------------------------------------------

      y_lo = 0.01d0
      Gauss_maxT = 2000.0d0
      Gauss_T_width = 0.0012d0
      Gaussian_Spec_width = 0.001d0

!      write(6,*)" made it to initdata"
      if (bathID.lt.1 .or. bathID.gt.nspecies) then
         call bl_pd_abort()
      endif

      ! Set up hot air composition
      Yl(:,1) = 0.0d0
      Yl(bathID-1,1) = 0.767
      Yl(oxidID-1,1) = 0.233

      Patm = pamb / pphys_getP1atm_MKS()

      do j = lo(2), hi(2)
        y = (float(j)+.5d0)*delta(2)+domnlo(2)
        do i = lo(1), hi(1)
          x = (float(i)+.5d0)*delta(1)+domnlo(1)
               
          ! Start with mixture fraction tanh between 0 and 1
          tanhval = 0.5d0*(1.0d0+TANH((x-midtanh)/widthtanh))
          call set_Y_from_Ksi(tanhval,Yl(0:nspecies-1,2))

          ! Set background T
          scal(i,j,Temp) = T_bc(1)

          ! Setup hot air region in the mixing layer
          rad = SQRT((x-splitx)**2.0+(MIN(y,y_lo)-y_lo)**2.0)
          Gaussian_T = EXP(-(rad)**2.0/(2.0*(Gauss_T_width + Gauss_T_width * (y-y_lo)/0.03)**2.0))
          scal(i,j,Temp) = T_bc(1) + (Gauss_maxT-T_bc(1)) * Gaussian_T

          Gaussian_spec = EXP(-(rad)**2.0/(2.0*(Gaussian_Spec_width + Gaussian_Spec_width * (y-y_lo)/0.03)**2.0))
          do n = 0, nspecies-1
             scal(i,j,FirstSpec+n) = (1.0d0 - Gaussian_Spec) * Yl(n,2)  + Gaussian_spec * Yl(n,1)
          end do

          ! Maybe not necessary, but to be sure, normalized mass fraction ...
          sumspec = SUM(scal(i,j,FirstSpec:FirstSpec+nspecies-1))
          do n = 0, nspecies-1
             scal(i,j,FirstSpec+n) = scal(i,j,FirstSpec+n) / sumspec
          end do

          scal(i,j,Trac) = 0.d0

          vel(i,j,1) = 0.d0
          vel(i,j,2) = V_in 

        end do
      end do

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
