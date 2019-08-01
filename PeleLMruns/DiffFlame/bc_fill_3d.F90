#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <PeleLM_F.H>
#include <AMReX_ArrayLim.H>

module bc_fill_3d_module

  use amrex_fort_module, only : dim=>amrex_spacedim

  implicit none
  
  private
  
  public :: den_fill, adv_fill, &
            temp_fill, rhoh_fill, vel_fill, all_chem_fill, &
            xvel_fill, yvel_fill, zvel_fill, chem_fill, press_fill

contains

  subroutine den_fill (den,DIMS(den),domlo,domhi,delta, xlo,time,bc) bind(C, name="den_fill")
                       
    use network, only : nspecies
    use mod_Fvar_def, only : domnlo
    use user_defined_fcts_3d_module, only : bcfunction
              
    implicit none

    integer DIMDEC(den), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  den(DIMV(den))
      
    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:nspecies-1), T, h

    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(den)
    lo(2) = ARG_L2(den)
    lo(3) = ARG_L3(den)
    hi(1) = ARG_H1(den)
    hi(2) = ARG_H2(den)
    hi(3) = ARG_H3(den)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))
      
    call filcc (den,DIMS(den),domlo,domhi,delta,xlo,bc)

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,y,z,3,1,time,u,v,w,rho,Yl,T,h,delta,.false.)
                den(i,j,k) = rho
             enddo
          enddo
       enddo
    endif
      
  end subroutine den_fill

  subroutine adv_fill (adv,DIMS(adv),domlo,domhi,delta,xlo,time,bc) bind(C, name="adv_fill")

    implicit none

    integer    DIMDEC(adv)
    integer    domlo(dim), domhi(dim)
    REAL_T     delta(dim), xlo(dim), time
    REAL_T     adv(DIMV(adv))
    integer    bc(dim,2)

    integer    i,j,k
    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(adv)
    lo(2) = ARG_L2(adv)
    lo(3) = ARG_L3(adv)
    hi(1) = ARG_H1(adv)
    hi(2) = ARG_H2(adv)
    hi(3) = ARG_H3(adv)

    call filcc (adv,DIMS(adv),domlo,domhi,delta,xlo,bc)

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          do j = lo(2),hi(2)
             do i = lo(1), hi(1)
                adv(i,j,k) = 0.0d0
             enddo
          enddo
       enddo
    endif

  end subroutine adv_fill

  subroutine temp_fill  (temp,DIMS(temp),domlo,domhi,delta, xlo,time,bc) bind(C, name="temp_fill")

    use network, only : nspecies
    use mod_Fvar_def, only : domnlo
    use user_defined_fcts_3d_module, only : bcfunction
      
    implicit none

    integer DIMDEC(temp), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  temp(DIMV(temp))
  
    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:nspecies-1), T, h

    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(temp)
    lo(2) = ARG_L2(temp)
    lo(3) = ARG_L3(temp)
    hi(1) = ARG_H1(temp)
    hi(2) = ARG_H2(temp)
    hi(3) = ARG_H3(temp)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))
      
    call filcc (temp,DIMS(temp),domlo,domhi,delta,xlo,bc)

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,y,z,3,1,time,u,v,w,rho,Yl,T,h,delta,.false.)
                temp(i,j,k) = T
             enddo
          enddo
       enddo
    endif
      
  end subroutine temp_fill

  subroutine rhoh_fill  (rhoh,DIMS(rhoh),domlo,domhi,delta, xlo,time,bc) bind(C, name="rhoh_fill")

    use network, only : nspecies
    use mod_Fvar_def, only : domnlo
    use user_defined_fcts_3d_module, only : bcfunction
      
    implicit none

    integer DIMDEC(rhoh), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  rhoh(DIMV(rhoh))
      
    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:nspecies-1), T, h

    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(rhoh)
    lo(2) = ARG_L2(rhoh)
    lo(3) = ARG_L3(rhoh)
    hi(1) = ARG_H1(rhoh)
    hi(2) = ARG_H2(rhoh)
    hi(3) = ARG_H3(rhoh)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))
      
    call filcc (rhoh,DIMS(rhoh),domlo,domhi,delta,xlo,bc)

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,y,z,3,1,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoh(i,j,k) = rho*h
             enddo
          enddo
       enddo
    endif

  end subroutine rhoh_fill

  subroutine vel_fill  (vel,DIMS(vel),domlo,domhi,delta,xlo,time,bc) bind(C, name="vel_fill")

    use network, only : nspecies
    use mod_Fvar_def, only : domnlo
    use user_defined_fcts_3d_module, only : bcfunction
#if defined(BL_DO_FLCT)
    use turbinflow_module
    use probdata_module, only : do_flct
    use mod_Fvar_def, only : V_in
#endif    
    use probdata_module, only : splitx, xfrontw, turb_scale

    implicit none

    integer DIMDEC(vel), bc(dim,2,dim)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  vel(DIMV(vel),dim)
#if defined(BL_DO_FLCT)
    REAL_T  xx(vel_l1:vel_h1), yy(vel_l2:vel_h2), zz, vfluc(vel_l1:vel_h1,vel_l2:vel_h2,turb_ncomp)
#endif
    integer i, j, k
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:nspecies-1), T, h, eta1

    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(vel)
    lo(2) = ARG_L2(vel)
    lo(3) = ARG_L3(vel)
    hi(1) = ARG_H1(vel)
    hi(2) = ARG_H2(vel)
    hi(3) = ARG_H3(vel)

    call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),1), &
                DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,1))
    call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),2), &
                DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,2))
    call filcc (vel(ARG_L1(vel),ARG_L2(vel),ARG_L3(vel),3), &
                DIMS(vel),domlo,domhi,delta,xlo,bc(1,1,3))

    if ((lo(3).lt.domlo(3)) &
         .and. ( (bc(3,1,1).eq.EXT_DIR) &
         .or.    (bc(3,1,2).eq.EXT_DIR) &
         .or.    (bc(3,1,3).eq.EXT_DIR) ) ) then

       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)

                call bcfunction(x,y,z,3,1,time,u,v,w,rho,Yl,T,h,delta,.true.)

                if (bc(3,1,1).eq.EXT_DIR) vel(i,j,k,1) = u
                if (bc(3,1,2).eq.EXT_DIR) vel(i,j,k,2) = v
                if (bc(3,1,3).eq.EXT_DIR) vel(i,j,k,3) = w

             enddo
          enddo
       enddo
    
#if defined(BL_DO_FLCT)
       if (do_flct) then
          do i = vel_l1,vel_h1
             xx(i) = (float(i)+.5)*delta(1)+domnlo(1)            
          enddo
          do j = vel_l2,vel_h2
             yy(j) = (float(j)+.5)*delta(2)+domnlo(2)            
          enddo
          zz = time*V_in
          vfluc = 0.d0
          call get_turbstate(vel_l1,vel_l2,vel_h1,vel_h2,xx,yy,zz,vfluc)
          do k = vel_l3, domlo(3)-1
             do j = vel_l2,vel_h2
                eta1 = 0.5d0*(1.d0 - TANH(2.d0*(ABS(yy(j))-splitx)/xfrontw))
                if (bc(3,1,1).eq.EXT_DIR) vel(:,j,k,1) = vel(:,j,k,1) + eta1*vfluc(:,j,1)*turb_scale
                if (bc(3,1,2).eq.EXT_DIR) vel(:,j,k,2) = vel(:,j,k,2) + eta1*vfluc(:,j,2)*turb_scale
                if (bc(3,1,3).eq.EXT_DIR) vel(:,j,k,3) = vel(:,j,k,3) + eta1*vfluc(:,j,3)*turb_scale
             enddo
          enddo
       endif
#endif
    endif

  end subroutine vel_fill

  subroutine xvel_fill (xvel,DIMS(xvel),domlo,domhi,delta, &
       xlo,time,bc) bind(C, name="xvel_fill")

    implicit none

    integer DIMDEC(xvel), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  xvel(DIMV(xvel))

    call bl_pd_abort('should not be in xvel_fill')
      
  end subroutine xvel_fill

  subroutine yvel_fill (yvel,DIMS(yvel),domlo,domhi,delta,&
       xlo,time,bc)  bind(C, name="yvel_fill")

    implicit none

    integer DIMDEC(yvel), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  yvel(DIMV(yvel))

    call bl_pd_abort('should not be in yvel_fill')

  end subroutine yvel_fill

  subroutine zvel_fill (zvel,DIMS(zvel),domlo,domhi,delta, &
       xlo,time,bc) bind(C, name="zvel_fill")

    implicit none

    integer DIMDEC(zvel), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  zvel(DIMV(zvel))

    call bl_pd_abort('should not be in zvel_fill')

  end subroutine zvel_fill
  
  subroutine all_chem_fill(rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,time,bc) bind(C, name="all_chem_fill")

    use network,  only: nspecies
    use mod_Fvar_def, only : domnlo
    use user_defined_fcts_3d_module, only : bcfunction

    implicit none
      
    integer DIMDEC(rhoY), bc(dim,2)
    integer domlo(dim), domhi(dim)
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  rhoY(DIMV(rhoY),nspecies)

    integer i, j, k, n
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(nspecies), T, h

    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(rhoY)
    lo(2) = ARG_L2(rhoY)
    lo(3) = ARG_L3(rhoY)
    hi(1) = ARG_H1(rhoY)
    hi(2) = ARG_H2(rhoY)
    hi(3) = ARG_H3(rhoY)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))
      
    do n = 1,nspecies
       call filcc (rhoY(lo(1),lo(2),lo(3),n), &
                   DIMS(rhoY),domlo,domhi,delta,xlo,bc)
    end do

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,y,z,3,1,time,u,v,w,rho,Yl,T,h,delta,.false.)
                do n = 1,nspecies
                   rhoY(i,j,k,n) = rho*Yl(n)
                end do
             enddo
          enddo
       enddo
    endif
      
  end subroutine all_chem_fill

  subroutine chem_fill  (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,time,bc,id) bind(C, name="chem_fill")

    use network, only : nspecies
    use mod_Fvar_def, only : domnlo
    use user_defined_fcts_3d_module, only : bcfunction
      
    implicit none

    integer DIMDEC(rhoY), bc(dim,2)
    integer domlo(dim), domhi(dim), id
    REAL_T  delta(dim), xlo(dim), time
    REAL_T  rhoY(DIMV(rhoY))
      
    integer i, j, k
    integer ilo, ihi, jlo, jhi, klo, khi
    REAL_T  z, y, x
    REAL_T  u, v, w, rho, Yl(0:nspecies-1), T, h

    integer lo(dim), hi(dim)

    lo(1) = ARG_L1(rhoY)
    lo(2) = ARG_L2(rhoY)
    lo(3) = ARG_L3(rhoY)
    hi(1) = ARG_H1(rhoY)
    hi(2) = ARG_H2(rhoY)
    hi(3) = ARG_H3(rhoY)

    ilo = max(lo(1),domlo(1))
    jlo = max(lo(2),domlo(2))
    klo = max(lo(3),domlo(3))
    ihi = min(hi(1),domhi(1))
    jhi = min(hi(2),domhi(2))
    khi = min(hi(3),domhi(3))
      
    call filcc (rhoY,DIMS(rhoY),domlo,domhi,delta,xlo,bc)

    if (bc(3,1).eq.EXT_DIR.and.lo(3).lt.domlo(3)) then
       do k = lo(3), domlo(3)-1
          z = (float(k)+.5)*delta(3)+domnlo(3)
          do j = lo(2),hi(2)
             y = (float(j)+.5)*delta(2)+domnlo(2)
             do i = lo(1), hi(1)
                x = (float(i)+.5)*delta(1)+domnlo(1)
                call bcfunction(x,y,z,3,1,time,u,v,w,rho,Yl,T,h,delta,.false.)
                rhoY(i,j,k) = rho*Yl(id)
             enddo
          enddo
       enddo
    endif
      
  end subroutine chem_fill

  subroutine press_fill  (p,DIMS(p),domlo,domhi,dx,xlo,time,bc) bind(C, name="press_fill")

    implicit none

    integer    DIMDEC(p)
    integer    domlo(dim), domhi(dim)
    REAL_T     dx(dim), xlo(dim), time
    REAL_T     p(DIMV(p))
    integer    bc(dim,2)

    integer    i, j, k
    integer    ilo, ihi, jlo, jhi, klo, khi
    logical    fix_xlo, fix_xhi, fix_ylo, fix_yhi, fix_zlo, fix_zhi
    logical    per_xlo, per_xhi, per_ylo, per_yhi, per_zlo, per_zhi

    fix_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .ne. INT_DIR)
    per_xlo = (ARG_L1(p) .lt. domlo(1)) .and. (bc(1,1) .eq. INT_DIR)
    fix_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .ne. INT_DIR)
    per_xhi = (ARG_H1(p) .gt. domhi(1)) .and. (bc(1,2) .eq. INT_DIR)
    fix_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .ne. INT_DIR)
    per_ylo = (ARG_L2(p) .lt. domlo(2)) .and. (bc(2,1) .eq. INT_DIR)
    fix_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .ne. INT_DIR)
    per_yhi = (ARG_H2(p) .gt. domhi(2)) .and. (bc(2,2) .eq. INT_DIR)
    fix_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .ne. INT_DIR)
    per_zlo = (ARG_L3(p) .lt. domlo(3)) .and. (bc(3,1) .eq. INT_DIR)
    fix_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .ne. INT_DIR)
    per_zhi = (ARG_H3(p) .gt. domhi(3)) .and. (bc(3,2) .eq. INT_DIR)

    ilo = max(ARG_L1(p),domlo(1))
    jlo = max(ARG_L2(p),domlo(2))
    klo = max(ARG_L3(p),domlo(3))
    ihi = min(ARG_H1(p),domhi(1))
    jhi = min(ARG_H2(p),domhi(2))
    khi = min(ARG_H3(p),domhi(3))

    !***************
    !  SETTING BL_XLO
    !***************

    if (fix_xlo) then
       do i = ARG_L1(p), domlo(1)-1
          do k = klo, khi
             do j = jlo,jhi
                p(i,j,k) = p(ilo,j,k)
             end do
          end do
       end do

       if (fix_ylo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = klo, khi
                   p(i,j,k) = p(ilo,jlo,k)
                end do
             end do
          end do

          if (fix_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jlo,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jlo,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jlo,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jlo,k)
                   end do
                end do
             end do
          end if
       end if

       if (fix_yhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = klo, khi
                   p(i,j,k) = p(ilo,jhi,k)
                end do
             end do
          end do
          if (fix_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jhi,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,jhi,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jhi,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,jhi,k)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo, jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,klo)
                end do
             end do
          end do
          if (per_ylo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,j,klo)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ilo,j,klo)
                   end do
                end do
             end do
          end if

       end if

       if (fix_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo, jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,khi)
                end do
             end do
          end do
          if (per_ylo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,j,khi)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ilo,j,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_ylo) then
          do i = ARG_L1(p), domlo(1)-1
             do k = klo,khi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do i = ARG_L1(p), domlo(1)-1
             do k = klo,khi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo,jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = jlo,jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_ylo .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_ylo .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ilo,j,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_XHI
    !*****************************************************************************

    if (fix_xhi) then
       do i = domhi(1)+1, ARG_H1(p)
          do k = klo, khi
             do j = jlo,jhi
                p(i,j,k) = p(ihi,j,k)
             end do
          end do
       end do

       if (fix_ylo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = klo, khi
                   p(i,j,k) = p(ihi,jlo,k)
                end do
             end do
          end do

          if (fix_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jlo,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jlo,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jlo,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jlo,k)
                   end do
                end do
             end do
          end if
       end if
       if (fix_yhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = klo, khi
                   p(i,j,k) = p(ihi,jhi,k)
                end do
             end do
          end do
          if (fix_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jhi,klo)
                   end do
                end do
             end do
          else if (per_zlo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,jhi,k)
                   end do
                end do
             end do
          end if
          if (fix_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jhi,khi)
                   end do
                end do
             end do
          else if (per_zhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,jhi,k)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo, jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,klo)
                end do
             end do
          end do
          if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,j,klo)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(ihi,j,klo)
                   end do
                end do
             end do
          end if

       end if

       if (fix_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo, jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,khi)
                end do
             end do
          end do
          if (per_ylo) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,j,khi)
                   end do
                end do
             end do
          end if
          if (per_yhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(ihi,j,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_ylo) then
          do i = domhi(1)+1, ARG_H1(p)
             do k = klo,khi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do k = klo,khi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo,jhi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = jlo,jhi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if


       if (per_ylo .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_ylo .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

       if (per_yhi .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(ihi,j,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_YLO
    !*****************************************************************************

    if (fix_ylo) then
       do j = ARG_L2(p), domlo(2)-1
          do k = klo, khi
             do i = ilo, ihi
                p(i,j,k) = p(i,jlo,k)
             end do
          end do
       end do

       if (fix_zlo) then
          do j = ARG_L2(p), domlo(2)-1
             do k = ARG_L3(p), domlo(3)-1
                do i = ilo, ihi
                   p(i,j,k) = p(i,jlo,klo)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jlo,klo)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jlo,klo)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zhi) then
          do j = ARG_L2(p), domlo(2)-1
             do k = domhi(3)+1, ARG_H3(p)
                do i = ilo, ihi
                   p(i,j,k) = p(i,jlo,khi)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jlo,khi)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jlo,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_xlo) then
          do j = ARG_L2(p), domlo(2)-1
             do k = klo,khi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do j = ARG_L2(p), domlo(2)-1
             do k = klo,khi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do j = ARG_L2(p), domlo(2)-1
             do i = ilo,ihi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do j = ARG_L2(p), domlo(2)-1
             do i = ilo,ihi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if


       if (per_xlo .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = ARG_L2(p), domlo(2)-1
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jlo,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_YHI
    !*****************************************************************************

    if (fix_yhi) then
       do j = domhi(2)+1, ARG_H2(p)
          do k = klo, khi
             do i = ilo, ihi
                p(i,j,k) = p(i,jhi,k)
             end do
          end do
       end do

       if (fix_zlo) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = ARG_L3(p), domlo(3)-1
                do i = ilo, ihi
                   p(i,j,k) = p(i,jhi,klo)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jhi,klo)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = ARG_L3(p), domlo(3)-1
                      p(i,j,k) = p(i,jhi,klo)
                   end do
                end do
             end do
          end if
       end if

       if (fix_zhi) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = domhi(3)+1, ARG_H3(p)
                do i = ilo, ihi
                   p(i,j,k) = p(i,jhi,khi)
                end do
             end do
          end do
          if (per_xlo) then
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jhi,khi)
                   end do
                end do
             end do
          end if
          if (per_xhi) then
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   do k = domhi(3)+1, ARG_H3(p)
                      p(i,j,k) = p(i,jhi,khi)
                   end do
                end do
             end do
          end if
       end if

       if (per_xlo) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = klo,khi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do j = domhi(2)+1, ARG_H2(p)
             do k = klo,khi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_zlo) then
          do j = domhi(2)+1, ARG_H2(p)
             do i = ilo,ihi
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if
       if (per_zhi) then
          do j = domhi(2)+1, ARG_H2(p)
             do i = ilo,ihi
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_zlo) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_zhi) then
          do i = ARG_L1(p), domlo(1)-1
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zlo) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = ARG_L3(p), domlo(3)-1
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_zhi) then
          do i = domhi(1)+1, ARG_H1(p)
             do j = domhi(2)+1, ARG_H2(p)
                do k = domhi(3)+1, ARG_H3(p)
                   p(i,j,k) = p(i,jhi,k)
                end do
             end do
          end do
       end if

    end if

    !*****************************************************************************
    ! SETTING BL_ZLO
    !*****************************************************************************

    if (fix_zlo) then
       do k = ARG_L3(p), domlo(3)-1
          do j = jlo, jhi
             do i = ilo, ihi
                p(i,j,k) = p(i,j,klo)
             end do
          end do
       end do

       if (per_xlo) then
          do k = ARG_L3(p), domlo(3)-1
             do j = jlo,jhi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do k = ARG_L3(p), domlo(3)-1
             do j = jlo,jhi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_ylo) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ilo,ihi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ilo,ihi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_ylo) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_yhi) then
          do k = ARG_L3(p), domlo(3)-1
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_ylo) then
          do k = ARG_L3(p), domlo(3)-1
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,klo)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_yhi) then
          do k = ARG_L3(p), domlo(3)-1
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
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
       do k = domhi(3)+1, ARG_H3(p)
          do j = jlo, jhi
             do i = ilo, ihi
                p(i,j,k) = p(i,j,khi)
             end do
          end do
       end do

       if (per_xlo) then
          do k = domhi(3)+1, ARG_H3(p)
             do j = jlo,jhi
                do i = ARG_L1(p), domlo(1)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if
       if (per_xhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do j = jlo,jhi
                do i = domhi(1)+1, ARG_H1(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_ylo) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ilo,ihi
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if
       if (per_yhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ilo,ihi
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if


       if (per_xlo .and. per_ylo) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ARG_L1(p), domlo(1)-1
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_xlo .and. per_yhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = ARG_L1(p), domlo(1)-1
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_ylo) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = domhi(1)+1, ARG_H1(p)
                do j = ARG_L2(p), domlo(2)-1
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

       if (per_xhi .and. per_yhi) then
          do k = domhi(3)+1, ARG_H3(p)
             do i = domhi(1)+1, ARG_H1(p)
                do j = domhi(2)+1, ARG_H2(p)
                   p(i,j,k) = p(i,j,khi)
                end do
             end do
          end do
       end if

    end if

  end subroutine press_fill

end module bc_fill_3d_module
