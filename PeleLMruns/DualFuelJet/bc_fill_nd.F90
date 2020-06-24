#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>
#include "mechanism.h"

module bc_fill_nd_module

  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_filcc_module, only : amrex_filccn
  use mod_Fvar_def, only : domnlo
  use user_defined_fcts_nd_module, only : bcfunction

  implicit none
  
  private
  
  public :: den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            xvel_fill, yvel_fill, zvel_fill, chem_fill, press_fill

  integer, parameter :: nspecies = NUM_SPECIES
  
contains

  subroutine den_fill (den, d_lo, d_hi, &
                       domlo, domhi, delta, &
                       xlo, time, bc)&
                       bind(C, name="den_fill")
      
    implicit none

! In/Out      
    integer :: d_lo(3), d_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3)) :: den

    integer i, j, k
    REAL_T  x(3)
    REAL_T  vel(3), rho, Yl(0:nspecies-1), T, h

    call amrex_filccn (d_lo, d_hi, den, d_lo, d_hi, 1, domlo, domhi, delta, xlo, bc)

    if (bc(3,1).eq.EXT_DIR.and.d_lo(3).lt.domlo(3)) then
       do k = d_lo(3), domlo(3)-1
          x(3) = (float(k)+.5)*delta(3)+domnlo(3)
          do j = d_lo(2),d_hi(2)
             x(2) = (float(j)+.5)*delta(2)+domnlo(2)
             do i = d_lo(1), d_hi(1)
                x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,delta,3,1,time,.false.,vel,rho,Yl,T,h)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif
      
  end subroutine den_fill

  subroutine adv_fill (adv, a_lo, a_hi, &
                       domlo, domhi, delta, &
                       xlo, time, bc)&
                       bind(C, name="adv_fill")
      
    implicit none

! In/Out      
    integer :: a_lo(3), a_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(a_lo(1):a_hi(1),a_lo(2):a_hi(2),a_lo(3):a_hi(3)) :: adv

    integer    i,j,k

    call amrex_filccn (a_lo, a_hi, adv, a_lo, a_hi, 1, domlo, domhi, delta, xlo, bc)

    if (bc(3,1).eq.EXT_DIR.and.a_lo(3).lt.domlo(3)) then
       do k = a_lo(3), domlo(3)-1
          do j = a_lo(2),a_hi(2)
             do i = a_lo(1), a_hi(1)
                adv(i,j,k) = 0.0d0
             enddo
          enddo
       enddo
    endif

  end subroutine adv_fill

  subroutine temp_fill (temp, t_lo, t_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="temp_fill")
      
    implicit none

! In/Out      
    integer :: t_lo(3), t_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) :: temp

    integer i, j, k
    REAL_T  x(3)
    REAL_T  vel(3), rho, Yl(0:nspecies-1), T, h

    call amrex_filccn (t_lo, t_hi, temp, t_lo, t_hi, 1, domlo, domhi, delta, xlo, bc)

    if (bc(3,1).eq.EXT_DIR.and.t_lo(3).lt.domlo(3)) then
       do k = t_lo(3), domlo(3)-1
          x(3) = (float(k)+.5)*delta(3)+domnlo(3)
          do j = t_lo(2),t_hi(2)
             x(2) = (float(j)+.5)*delta(2)+domnlo(2)
             do i = t_lo(1), t_hi(1)
                x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,delta,3,1,time,.false.,vel,rho,Yl,T,h)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif
      
  end subroutine temp_fill

  subroutine rhoh_fill (rhoh, r_lo, r_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="rhoh_fill")
      
    implicit none

! In/Out      
    integer :: r_lo(3), r_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: rhoh

    integer i, j, k
    REAL_T  x(3)
    REAL_T  vel(3), rho, Yl(0:nspecies-1), T, h

    call amrex_filccn (r_lo, r_hi, rhoh, r_lo, r_hi, 1, domlo, domhi, delta, xlo, bc)

    if (bc(3,1).eq.EXT_DIR.and.r_lo(3).lt.domlo(3)) then
       do k = r_lo(3), domlo(3)-1
          x(3) = (float(k)+.5)*delta(3)+domnlo(3)
          do j = r_lo(2),r_hi(2)
             x(2) = (float(j)+.5)*delta(2)+domnlo(2)
             do i = r_lo(1), r_hi(1)
                x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,delta,3,1,time,.false.,vel,rho,Yl,T,h)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

  end subroutine rhoh_fill

  subroutine vel_fill (vel, v_lo, v_hi, &
                       domlo, domhi, delta, &
                       xlo, time, bc)&
                       bind(C, name="vel_fill")
      
#if defined(BL_DO_FLCT)
    use turbinflow_module
    use probdata_module, only : do_flct
    use mod_Fvar_def, only : V_in
#endif    
    use probdata_module, only : splitx, xfrontw, turb_scale

    implicit none

! In/Out      
    integer :: v_lo(3), v_hi(3)
    integer :: bc(dim,2,3)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3), dim) :: vel

#if defined(BL_DO_FLCT)
    REAL_T  xx(v_lo(1):v_hi(1)), yy(v_lo(2):v_hi(2)), zz, vfluc(v_lo(1):v_hi(1),v_lo(2):v_hi(2),turb_ncomp)
#endif
    integer i, j, k
    REAL_T  x(3)
    REAL_T  u(3), rho, Yl(0:nspecies-1), T, h, eta1

    call amrex_filccn (v_lo, v_hi, vel(v_lo(1),v_lo(2),v_lo(3),1), v_lo, v_hi, 1, domlo, domhi, delta, xlo, bc(1,1,1))
    call amrex_filccn (v_lo, v_hi, vel(v_lo(1),v_lo(2),v_lo(3),2), v_lo, v_hi, 1, domlo, domhi, delta, xlo, bc(1,1,2))
#if BL_SPACEDIM > 2
    call amrex_filccn (v_lo, v_hi, vel(v_lo(1),v_lo(2),v_lo(3),3), v_lo, v_hi, 1, domlo, domhi, delta, xlo, bc(1,1,3))
#endif

    if ((v_lo(3).lt.domlo(3)) &
         .and. ( (bc(3,1,1).eq.EXT_DIR) &
         .or.    (bc(3,1,2).eq.EXT_DIR) &
         .or.    (bc(3,1,3).eq.EXT_DIR) ) ) then

       do k = v_lo(3), domlo(3)-1
          x(3) = (float(k)+.5)*delta(3)+domnlo(3)
          do j = v_lo(2),v_hi(2)
             x(2) = (float(j)+.5)*delta(2)+domnlo(2)
             do i = v_lo(1), v_hi(1)
                x(1) = (float(i)+.5)*delta(1)+domnlo(1)

                call bcfunction(x,delta,3,1,time,.true.,u,rho,Yl,T,h)

                if (bc(3,1,1).eq.EXT_DIR) vel(i,j,k,1) = u(1)
                if (bc(3,1,2).eq.EXT_DIR) vel(i,j,k,2) = u(2)
                if (bc(3,1,3).eq.EXT_DIR) vel(i,j,k,3) = u(3)

             enddo
          enddo
       enddo
    
#if defined(BL_DO_FLCT)
       if (do_flct) then
          do i = v_lo(1),v_hi(1)
             xx(i) = (float(i)+.5)*delta(1)+domnlo(1)            
          enddo
          do j = v_lo(2),v_hi(2)
             yy(j) = (float(j)+.5)*delta(2)+domnlo(2)            
          enddo
          zz = time*V_in
          vfluc = 0.d0
          call get_turbstate(v_lo(1),v_lo(2),v_hi(1),v_hi(2),xx,yy,zz,vfluc)
          do k = v_lo(3), domlo(3)-1
             do j = v_lo(2),v_hi(2)
              do i = v_lo(1),v_hi(1)
                eta1 = 0.5d0*(1.d0 - TANH(2.d0*(SQRT(yy(j)**2+xx(i)**2)-splitx)/xfrontw))
                if (bc(3,1,1).eq.EXT_DIR) vel(i,j,k,1) = vel(i,j,k,1) + eta1*vfluc(i,j,1)*turb_scale
                if (bc(3,1,2).eq.EXT_DIR) vel(i,j,k,2) = vel(i,j,k,2) + eta1*vfluc(i,j,2)*turb_scale
                if (bc(3,1,3).eq.EXT_DIR) vel(i,j,k,3) = vel(i,j,k,3) + eta1*vfluc(i,j,3)*turb_scale

               enddo
             enddo
          enddo
       endif
#endif
    endif

  end subroutine vel_fill

  subroutine xvel_fill (vel, v_lo, v_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="xvel_fill")

    implicit none

! In/Out      
    integer :: v_lo(3), v_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vel

    call bl_pd_abort('should not be in xvel_fill')
      
  end subroutine xvel_fill

  subroutine yvel_fill (vel, v_lo, v_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="yvel_fill")
 
    implicit none

! In/Out      
    integer :: v_lo(3), v_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vel

    call bl_pd_abort('should not be in yvel_fill')

  end subroutine yvel_fill

  subroutine zvel_fill (vel, v_lo, v_hi, &
                        domlo, domhi, delta, &
                        xlo, time, bc)&
                        bind(C, name="zvel_fill")
      
    implicit none

! In/Out      
    integer :: v_lo(3), v_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(v_lo(1):v_hi(1),v_lo(2):v_hi(2),v_lo(3):v_hi(3)) :: vel

    call bl_pd_abort('should not be in zvel_fill')

  end subroutine zvel_fill

  subroutine all_chem_fill (rhoY, r_lo, r_hi, &
                            domlo, domhi, delta, &
                            xlo, time, bc)&
                            bind(C, name="all_chem_fill")

    implicit none

! In/Out      
    integer :: r_lo(3), r_hi(3)
    integer :: bc(dim,2)
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),NUM_SPECIES) :: rhoY

    integer i, j, k, n
    REAL_T  x(3)
    REAL_T  vel(3), rho, Yl(nspecies), T, h

    do n = 1,nspecies
       call filcc (r_lo, r_hi, rhoY(r_lo(1),r_lo(2),r_lo(3),n), &
                   r_lo, r_hi, 1, domlo, domhi, delta, xlo, bc)
    end do

    if (bc(3,1).eq.EXT_DIR.and.r_lo(3).lt.domlo(3)) then
       do k = r_lo(3), domlo(3)-1
          x(3) = (float(k)+.5)*delta(3)+domnlo(3)
          do j = r_lo(2),r_hi(2)
             x(2) = (float(j)+.5)*delta(2)+domnlo(2)
             do i = r_lo(1), r_hi(1)
                x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,delta,3,1,time,.false.,vel,rho,Yl,T,h)
                do n = 1,nspecies
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
             enddo
          enddo
       enddo
    endif
      
  end subroutine all_chem_fill

  subroutine chem_fill (rhoY, r_lo, r_hi, &
                            domlo, domhi, delta, &
                            xlo, time, bc, id)&
                            bind(C, name="chem_fill")

    implicit none

! In/Out      
    integer :: r_lo(3), r_hi(3)
    integer :: bc(dim,2), id
    integer :: domlo(3), domhi(3)
    REAL_T  :: delta(3), xlo(3), time
    REAL_T, dimension(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3)) :: rhoY

    integer i, j, k
    REAL_T  x(3)
    REAL_T  vel(3), rho, Yl(0:nspecies-1), T, h

    call amrex_filccn (r_lo, r_hi, rhoY, r_lo, r_hi, 1, domlo, domhi, delta, xlo, bc)

    if (bc(3,1).eq.EXT_DIR.and.r_lo(3).lt.domlo(3)) then
       do k = r_lo(3), domlo(3)-1
          x(3) = (float(k)+.5)*delta(3)+domnlo(3)
          do j = r_lo(2),r_hi(2)
             x(2) = (float(j)+.5)*delta(2)+domnlo(2)
             do i = r_lo(1), r_hi(1)
                x(1) = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,delta,3,1,time,.false.,vel,rho,Yl,T,h)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif
      
  end subroutine chem_fill

   subroutine press_fill (p, p_lo, p_hi, &
                          domlo, domhi, delta, &
                          xlo, time, bc)&
                          bind(C, name="press_fill")

      implicit none

! In/Out
      integer :: p_lo(3), p_hi(3)
      integer :: bc(dim,2)
      integer :: domlo(3), domhi(3)
      REAL_T  :: delta(3), xlo(3), time
      REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)) :: p

! Local
      integer :: i, j, k
      integer :: ilo, ihi, jlo, jhi, klo, khi
      logical :: fix_xlo, fix_xhi, fix_ylo, fix_yhi, fix_zlo, fix_zhi
      logical :: per_xlo, per_xhi, per_ylo, per_yhi, per_zlo, per_zhi

      fix_xlo = (p_lo(1) < domlo(1)) .and. (bc(1,1) /= INT_DIR)
      per_xlo = (p_lo(1) < domlo(1)) .and. (bc(1,1) == INT_DIR)
      fix_xhi = (p_hi(1) > domhi(1)) .and. (bc(1,2) /= INT_DIR)
      per_xhi = (p_hi(1) > domhi(1)) .and. (bc(1,2) == INT_DIR)
#if ( AMREX_SPACEDIM >= 2 )      
      fix_ylo = (p_lo(2) < domlo(2)) .and. (bc(2,1) /= INT_DIR)
      per_ylo = (p_lo(2) < domlo(2)) .and. (bc(2,1) == INT_DIR)
      fix_yhi = (p_hi(2) > domhi(2)) .and. (bc(2,2) /= INT_DIR)
      per_yhi = (p_hi(2) > domhi(2)) .and. (bc(2,2) == INT_DIR)
#if ( AMREX_SPACEDIM == 3 )      
      fix_zlo = (p_lo(3) < domlo(3)) .and. (bc(3,1) /= INT_DIR)
      per_zlo = (p_lo(3) < domlo(3)) .and. (bc(3,1) == INT_DIR)
      fix_zhi = (p_hi(3) > domhi(3)) .and. (bc(3,2) /= INT_DIR)
      per_zhi = (p_hi(3) > domhi(3)) .and. (bc(3,2) == INT_DIR)
#endif
#endif

      ilo = max(p_lo(1),domlo(1))
      ihi = min(p_hi(1),domhi(1))
      jlo = p_lo(2)
      klo = p_lo(3)
      jhi = p_hi(2)
      khi = p_hi(3)
#if ( AMREX_SPACEDIM >= 2 )      
      jlo = max(p_lo(2),domlo(2))
      jhi = min(p_hi(2),domhi(2))
#if ( AMREX_SPACEDIM == 3 )      
      klo = max(p_lo(3),domlo(3))
      khi = min(p_hi(3),domhi(3))
#endif
#endif

!***************
!  SETTING BL_XLO
!***************

      if (fix_xlo) then
         do i = p_lo(1), domlo(1)-1
            do k = p_lo(3), p_hi(3)
               do j = p_lo(2), p_hi(2)
                  p(i,j,k) = p(ilo,j,k)
               end do 
            end do
         end do

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_ylo) then
            do i = p_lo(1), domlo(1)-1
                 do j = p_lo(2), domlo(2)-1
                    do k = klo, khi
                       p(i,j,k) = p(ilo,jlo,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jlo,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jlo,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jlo,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jlo,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_yhi) then
            do i = p_lo(1), domlo(1)-1
                 do j = domhi(2)+1, p_hi(2)
                    do k = klo, khi
                       p(i,j,k) = p(ilo,jhi,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jhi,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,jhi,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jhi,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,jhi,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (fix_zlo) then
            do i = p_lo(1), domlo(1)-1
                 do j = jlo, jhi
                    do k = p_lo(3), domlo(3)-1
                       p(i,j,k) = p(ilo,j,klo)
                    end do
                 end do
            end do
              if (per_ylo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,j,klo)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ilo,j,klo)
                       end do
                    end do
                 end do
              end if

         end if

         if (fix_zhi) then
            do i = p_lo(1), domlo(1)-1
                 do j = jlo, jhi
                    do k = domhi(3)+1, p_hi(3)
                       p(i,j,k) = p(ilo,j,khi)
                    end do
                 end do
            end do
              if (per_ylo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,j,khi)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
                 do i = p_lo(1), domlo(1)-1
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ilo,j,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif
 
#if ( AMREX_SPACEDIM >= 2 )      
         if (per_ylo) then
               do i = p_lo(1), domlo(1)-1
                  do k = klo,khi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do i = p_lo(1), domlo(1)-1
                  do k = klo,khi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
#endif
 
#if ( AMREX_SPACEDIM == 3 )      
         if (per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = jlo,jhi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = jlo,jhi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ilo,j,k)
                     end do
                  end do
               end do
         end if
#endif

      end if     ! End if on fix_xlo

!*****************************************************************************
! SETTING BL_XHI
!*****************************************************************************

      if (fix_xhi) then
         do i = domhi(1)+1, p_hi(1)
            do k = klo, khi
               do j = jlo,jhi
                  p(i,j,k) = p(ihi,j,k)
               end do
            end do
         end do

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_ylo) then
            do i = domhi(1)+1, p_hi(1)
                 do j = p_lo(2), domlo(2)-1
                    do k = klo, khi
                       p(i,j,k) = p(ihi,jlo,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jlo,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jlo,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jlo,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jlo,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM >= 2 )      
         if (fix_yhi) then
            do i = domhi(1)+1, p_hi(1)
                 do j = domhi(2)+1, p_hi(2)
                    do k = klo, khi
                       p(i,j,k) = p(ihi,jhi,k)
                    end do
                 end do
            end do

#if ( AMREX_SPACEDIM == 3 )      
            if (fix_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jhi,klo)
                       end do
                    end do
                 end do
            else if (per_zlo) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,jhi,k)
                       end do
                    end do
                 end do
            end if
            if (fix_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jhi,khi)
                       end do
                    end do
                 end do
            else if (per_zhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,jhi,k)
                       end do
                    end do
                 end do
            end if
#endif
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (fix_zlo) then
            do i = domhi(1)+1, p_hi(1)
                 do j = jlo, jhi
                    do k = p_lo(3), domlo(3)-1
                       p(i,j,k) = p(ihi,j,klo)
                    end do
                 end do
            end do
              if (per_ylo) then
               do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,j,klo)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
               do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(ihi,j,klo)
                       end do
                    end do
                 end do
              end if
         end if

         if (fix_zhi) then
            do i = domhi(1)+1, p_hi(1)
                 do j = jlo, jhi
                    do k = domhi(3)+1, p_hi(3)
                       p(i,j,k) = p(ihi,j,khi)
                    end do
                 end do
            end do
              if (per_ylo) then
               do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,j,khi)
                       end do
                    end do
                 end do
              end if
              if (per_yhi) then
               do i = domhi(1)+1, p_hi(1)
                    do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(ihi,j,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif

#if ( AMREX_SPACEDIM >= 2 )      
         if (per_ylo) then
             do i = domhi(1)+1, p_hi(1)
                  do k = klo,khi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
             do i = domhi(1)+1, p_hi(1)
                  do k = klo,khi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (per_zlo) then
             do i = domhi(1)+1, p_hi(1)
                  do j = jlo,jhi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
              do i = domhi(1)+1, p_hi(1)
                  do j = jlo,jhi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
#endif

#if ( AMREX_SPACEDIM == 3 )      
         if (per_ylo .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_ylo .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if

         if (per_yhi .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(ihi,j,k)
                     end do
                  end do
               end do
         end if
#endif

      end if   ! End if on fix_xhi

!*****************************************************************************
! SETTING BL_YLO
!*****************************************************************************

#if ( AMREX_SPACEDIM >= 2 )
      if (fix_ylo) then
         do j = p_lo(2), domlo(2)-1
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jlo,k)
               end do
            end do
         end do

#if ( AMREX_SPACEDIM == 3 )
         if (fix_zlo) then
            do j = p_lo(2), domlo(2)-1
                 do k = p_lo(3), domlo(3)-1
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jlo,klo)
                    end do
                 end do
            end do
            if (per_xlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
            if (per_xhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,klo)
                     end do
                  end do
               end do
            end if
         end if

         if (fix_zhi) then
            do j = p_lo(2), domlo(2)-1
                 do k = domhi(3)+1, p_hi(3)
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jlo,khi)
                    end do
                 end do
            end do
              if (per_xlo) then
                 do i = p_lo(1), domlo(1)-1
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jlo,khi)
                       end do
                    end do
                 end do
              end if
              if (per_xhi) then
                 do i = domhi(1)+1, p_hi(1)
                    do j = p_lo(2), domlo(2)-1
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jlo,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif

         if (per_xlo) then
               do j = p_lo(2), domlo(2)-1
                  do k = klo,khi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = p_lo(2), domlo(2)-1
                  do k = klo,khi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

#if ( AMREX_SPACEDIM == 3 )
         if (per_zlo) then
               do j = p_lo(2), domlo(2)-1
                  do i = ilo,ihi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = p_lo(2), domlo(2)-1
                  do i = ilo,ihi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                  do j = p_lo(2), domlo(2)-1
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jlo,k)
                     end do
                  end do
               end do
         end if
#endif

      end if      ! End if on ylo
 
!*****************************************************************************
! SETTING BL_YHI
!*****************************************************************************

      if (fix_yhi) then
         do j = domhi(2)+1, p_hi(2)
            do k = klo, khi
               do i = ilo, ihi
                  p(i,j,k) = p(i,jhi,k)
               end do
            end do
         end do

#if ( AMREX_SPACEDIM == 3 )
         if (fix_zlo) then
            do j = domhi(2)+1, p_hi(2)
                 do k = p_lo(3), domlo(3)-1
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jhi,klo)
                    end do
                 end do
            end do
              if (per_xlo) then
                 do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(i,jhi,klo)
                       end do
                    end do
                 end do
              end if
              if (per_xhi) then
                 do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                       do k = p_lo(3), domlo(3)-1
                          p(i,j,k) = p(i,jhi,klo)
                       end do
                    end do
                 end do
              end if
         end if

         if (fix_zhi) then
            do j = domhi(2)+1, p_hi(2)
                 do k = domhi(3)+1, p_hi(3)
                    do i = ilo, ihi
                       p(i,j,k) = p(i,jhi,khi)
                    end do
                 end do
            end do
              if (per_xlo) then
                 do i = p_lo(1), domlo(1)-1
                  do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jhi,khi)
                       end do
                    end do
                 end do
              end if
              if (per_xhi) then
                 do i = domhi(1)+1, p_hi(1)
                  do j = domhi(2)+1, p_hi(2)
                       do k = domhi(3)+1, p_hi(3)
                          p(i,j,k) = p(i,jhi,khi)
                       end do
                    end do
                 end do
              end if
         end if
#endif

         if (per_xlo) then
               do j = domhi(2)+1, p_hi(2)
                  do k = klo,khi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do j = domhi(2)+1, p_hi(2)
                  do k = klo,khi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

#if ( AMREX_SPACEDIM == 3 )
         if (per_zlo) then
               do j = domhi(2)+1, p_hi(2)
                  do i = ilo,ihi
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
         if (per_zhi) then
               do j = domhi(2)+1, p_hi(2)
                  do i = ilo,ihi
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zlo) then
               do i = p_lo(1), domlo(1)-1
                do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_zhi) then
               do i = p_lo(1), domlo(1)-1
                do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zlo) then
               do i = domhi(1)+1, p_hi(1)
                do j = domhi(2)+1, p_hi(2)
                     do k = p_lo(3), domlo(3)-1
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_zhi) then
               do i = domhi(1)+1, p_hi(1)
                do j = domhi(2)+1, p_hi(2)
                     do k = domhi(3)+1, p_hi(3)
                        p(i,j,k) = p(i,jhi,k)
                     end do
                  end do
               end do
         end if
#endif

      end if     ! End if on yhi

!*****************************************************************************
! SETTING BL_ZLO
!*****************************************************************************

#if ( AMREX_SPACEDIM == 3 )
      if (fix_zlo) then
         do k = p_lo(3), domlo(3)-1
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,klo)
               end do
            end do
         end do

         if (per_xlo) then
               do k = p_lo(3), domlo(3)-1
                  do j = jlo,jhi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = p_lo(3), domlo(3)-1
                  do j = jlo,jhi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = p_lo(3), domlo(3)-1
                  do i = ilo,ihi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = p_lo(3), domlo(3)-1
                  do i = ilo,ihi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_ylo) then
               do k = p_lo(3), domlo(3)-1
                  do i = p_lo(1), domlo(1)-1
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = p_lo(3), domlo(3)-1
                  do i = p_lo(1), domlo(1)-1
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = p_lo(3), domlo(3)-1
                  do i = domhi(1)+1, p_hi(1)
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = p_lo(3), domlo(3)-1
                  do i = domhi(1)+1, p_hi(1)
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,klo)
                     end do
                  end do
               end do
         end if

      end if            

!*****************************************************************************
! SETTING BL_ZHI
!*****************************************************************************

      if (fix_zhi) then
         do k = domhi(3)+1, p_hi(3)
            do j = jlo, jhi
               do i = ilo, ihi
                  p(i,j,k) = p(i,j,khi)
               end do
            end do
       end do

         if (per_xlo) then
               do k = domhi(3)+1, p_hi(3)
                  do j = jlo,jhi
                     do i = p_lo(1), domlo(1)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_xhi) then
               do k = domhi(3)+1, p_hi(3)
                  do j = jlo,jhi
                     do i = domhi(1)+1, p_hi(1)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_ylo) then
               do k = domhi(3)+1, p_hi(3)
                  do i = ilo,ihi
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if
         if (per_yhi) then
               do k = domhi(3)+1, p_hi(3)
                  do i = ilo,ihi
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if


         if (per_xlo .and. per_ylo) then
               do k = domhi(3)+1, p_hi(3)
                  do i = p_lo(1), domlo(1)-1
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xlo .and. per_yhi) then
               do k = domhi(3)+1, p_hi(3)
                  do i = p_lo(1), domlo(1)-1
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_ylo) then
               do k = domhi(3)+1, p_hi(3)
                  do i = domhi(1)+1, p_hi(1)
                     do j = p_lo(2), domlo(2)-1
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

         if (per_xhi .and. per_yhi) then
               do k = domhi(3)+1, p_hi(3)
                  do i = domhi(1)+1, p_hi(1)
                     do j = domhi(2)+1, p_hi(2)
                        p(i,j,k) = p(i,j,khi)
                     end do
                  end do
               end do
         end if

      end if            
#endif
#endif

   end subroutine press_fill

end module bc_fill_nd_module
