#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module user_defined_fcts_2d_module

   use amrex_fort_module, only : dim=>amrex_spacedim
   use fuego_chemistry

implicit none

  private
  
  public :: bcfunction, zero_visc, set_Y_from_Phi, set_Y_from_ksi, set_Zst

contains
  


!-----------------------

  subroutine bcfunction(x,y,dir,norm,time,u,v,rho,Yl,T,h,dx,getuv) &
                        bind(C, name="bcfunction")

      use network,   only: nspecies
      use mod_Fvar_def, only : dv_control, tbase_control, V_in, f_flag_active_control, pamb
      use probdata_module, only : bcinit, rho_bc, Y_bc, T_bc, h_bc, v_bc, midtanh, widthtanh
      use PeleLM_F, only: pphys_getP1atm_MKS
      use PeleLM_2D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
      
      implicit none

      REAL_T x, y, time, u, v, rho, Yl(0:*), T, h, dx(dim)
      integer dir, norm  ! This specify the direction and orientation of the face
      logical getuv
      REAL_T :: tanhval, vbase, ksi
      REAL_T :: h_tmp(1), rho_tmp(1)

      integer n
      REAL_T :: Patm
      
      integer b(2)
      data  b / 1, 1 /

      if (.not. bcinit) then
         call bl_abort('Need to initialize boundary condition function')
      end if

      if ((dir == 2).and.(norm == 1)) then
        Patm = pamb / pphys_getP1atm_MKS()
        tanhval = 0.5d0*(1.0d0+TANH((x-midtanh)/widthtanh)) 
        ksi = tanhval
        call set_Y_from_Ksi(ksi,Yl(0:nspecies-1))
        T = T_bc(1)
        call pphys_RHOfromPTY(b, b, &
                             rho_tmp(1), DIMARG(b), DIMARG(b), &
                             T_bc(1),     DIMARG(b), DIMARG(b), &
                             Yl(0),   DIMARG(b), DIMARG(b), Patm)
        call pphys_HMIXfromTY(b, b, &
                             h_tmp(1), DIMARG(b), DIMARG(b), &
                             T_bc(1),   DIMARG(b), DIMARG(b), &
                             Yl(0), DIMARG(b), DIMARG(b))
        rho = rho_tmp(1)
        h = h_tmp(1)

        if (getuv .eqv. .TRUE.) then
            
          u = zero
          if (f_flag_active_control == 1 ) then               
            vbase =  V_in + (time-tbase_control)*dV_control
            v = vbase !+ tanhval * (vbase*rho_bc(1,2)/rho_bc(1,1) - vbase)
          else 
            v = V_in !+ tanhval * (V_in*rho_bc(1,2)/rho_bc(1,1) - V_in)
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
      use mod_Fvar_def, only : domnhi, domnlo
      
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

  subroutine set_Y_from_Phi(phi,Yt)bind(C, name="set_Y_from_Phi")
  
      use mod_Fvar_def, only : fuelID, oxidID, bathID
      use network,   only: nspecies
      use PeleLM_F,  only: pphys_get_spec_name2
      use probdata_module, only : H2_enrich

      implicit none

      REAL_T, INTENT(IN)  :: phi
      REAL_T, INTENT(OUT) :: Yt(nspecies)

      REAL_T :: a, a_hyd
      REAL_T :: Xt(nspecies)
      INTEGER ::  n
      CHARACTER(LEN=16) :: name

      Xt(:) = zero
      
!     Set "a" for computing X from phi
!     hc + a.O2 -> b.CO2 + c.H2O
      
!     ((1-alpha)*hc+alpha*H2) + a_hyd.O2 -> b_hyd.CO2 + c_hyd.H2O
!     a_hyd = (1-alpha)*a_hc + alpha*a_H2

      call pphys_get_spec_name2(name,fuelID)

      a = 0.d0
      if (name .eq. 'CH4') then
         a = 2.0d0
      else if (name .eq. 'C2H4') then
         a = 3.0d0
      else if (name .eq. 'C3H8') then
         a = 5.0d0
      else if (name .eq. 'C7H16') then
         a = 11.0d0
      else if (name .eq. 'H2') then
         a = 0.5d0
      else
         call bl_abort('setupbc: Unknown fuel type')
      end if

      !a_hyd = (1.0d0 - H2_enrich) * a + H2_enrich * 0.5d0
      a_hyd = a

      Xt(oxidID) = 1.d0/(1.d0 + phi/a_hyd  + 0.79d0/0.21d0)
      Xt(fuelID) = phi * Xt(oxidID) / a_hyd
      !Xt(fuelID) = phi * (1.0d0 - H2_enrich) * Xt(oxidID) / a_hyd
      !Xt(H2ID) = phi * H2_enrich * Xt(oxidID) / a_hyd
      !Xt(bathID) = 1.d0 - Xt(fuelID) - Xt(H2ID) - Xt(oxidID)
      Xt(bathID) = 1.d0 - Xt(fuelID) - Xt(oxidID)
      
      CALL CKXTY (Xt, Yt)
      
  end subroutine set_Y_from_Phi

  subroutine set_Y_from_Ksi(ksi,Yt)bind(C, name="set_Y_from_ksi")
  
      use network,   only: nspecies
      use probdata_module, only : Zst

      implicit none

      REAL_T, INTENT(IN)  :: ksi
      REAL_T, INTENT(OUT) :: Yt(nspecies)

      INTEGER ::  n
      REAL_T :: phi
      CHARACTER(LEN=16) :: name

      phi = ksi / ( 1.0 - min(ksi,0.999999) ) * ( 1.0 - Zst ) / Zst

      call set_Y_from_Phi(phi,Yt)
      
  end subroutine set_Y_from_ksi

  subroutine set_Zst( )bind(C, name="set_Zst")
  
      use network,   only: nspecies
      use PeleLM_F,  only: pphys_get_spec_name2
      use mod_Fvar_def, only : fuelID, oxidID
      use probdata_module, only : H2_enrich, Zst

      implicit none

      REAL_T :: a, a_hyd, mwt_fuel
      REAL_T :: MWT(nspecies)
      CHARACTER(LEN=16) :: name

      call pphys_get_spec_name2(name,fuelID)

      call CKWT(MWT)

      a = 0.d0
      if (name .eq. 'CH4') then
         a = 2.0d0
      else if (name .eq. 'C2H4') then
         a = 3.0d0
      else if (name .eq. 'C3H8') then
         a = 5.0d0
      else if (name .eq. 'C7H16') then
         a = 11.0d0
      else if (name .eq. 'H2') then
         a = 0.5d0
      else
         call bl_abort('setupbc: Unknown fuel type')
      end if

      !a_hyd = (1.0d0 - H2_enrich) * a + H2_enrich * 0.5d0
      !mwt_fuel = (1.0d0 - H2_enrich) * MWT(fuelID) + H2_enrich * MWT(H2ID)
      a_hyd = a
      mwt_fuel = MWT(fuelID)

      Zst = 1.0d0 / ( 1.0d0 + MWT(OxidID) / mwt_fuel * a_hyd / 0.233 ) 

      write(*,*) " Zst = ", Zst
      
  end subroutine set_Zst

end module user_defined_fcts_2d_module
