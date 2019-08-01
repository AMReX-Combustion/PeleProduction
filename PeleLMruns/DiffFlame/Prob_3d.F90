
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>
module prob_3D_module

  use amrex_fort_module, only : amrex_spacedim

  implicit none

  private

  public :: amrex_probinit, setupbc, init_data, set_Y_from_Phi

  
contains

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

    use network, only : nspecies, spec_names
    use PeleLM_F,  only : pphys_getP1atm_MKS
    use mod_Fvar_def, only : V_in

    use probdata_module, only : T_in, V_co, phi_in, T_co, &
                                splitx, xfrontw, fuel_N2_vol_percent, &
                                blobr, bloby, blobx, blobz, blobT, Tfrontw, turb_scale, &
                                domnlo, domnhi, adverr, diff_thresh, flametracval, H2_frac, &
                                hrr_thresh, max_diff_lev, max_hrr_lev, max_mix_lev, max_temp_lev, &
                                max_trac_lev, max_vort_lev, mix_thresh, OH_thresh, pertmag, &
                                pseudo_gravity, refine_nozzle, RO2_thresh, splity, standoff, &
                                stTh, T_switch, temperr, tempgrad, traceSpecVal, V_jet, &
                                vorterr, yfrontw, zmax_diff, zmax_mix, &
                                iH2, iO2, iN2, iCH4, iNC12H26

    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(amrex_spacedim), probhi(amrex_spacedim)

    integer i
    REAL_T pamb

    namelist /fortin/ vorterr, temperr, adverr, tempgrad, &
                      hrr_thresh, OH_thresh, RO2_thresh, &
                      flametracval, mix_thresh, zmax_mix,&
                      max_temp_lev, max_vort_lev, max_trac_lev, &
                      max_hrr_lev, max_mix_lev, &
                      traceSpecVal,phi_in,T_in, T_co, &
                      turb_scale, V_in, V_co, H2_frac,T_switch, &
                      standoff, pertmag, splitx, xfrontw, &
                      splity, yfrontw, blobx, bloby, blobz, blobr, &
                      blobT, Tfrontw, stTh, fuel_N2_vol_percent,V_jet,refine_nozzle, &
                      max_diff_lev, zmax_diff, diff_thresh
    namelist /heattransin/ pamb
#if defined(BL_DO_FLCT)
    namelist /flctin/ tstart_turb, forceInflow, numInflPlanesStore, forceLo, forceHi, &
         strmwse_dir, nCompInflow, flct_file
#endif

    !
    ! Build `probin' filename -- the name of file containing fortin namelist.
    !
    integer maxlen, isioproc
    parameter (maxlen=256)
    character probin*(maxlen)

    do i = 1,nspecies
       if (spec_names(i) .eq. "H2") iH2 = i
       if (spec_names(i) .eq. "O2") iO2 = i
       if (spec_names(i) .eq. "CH4") iCH4 = i
       if (spec_names(i) .eq. "N2") iN2 = i
       if (spec_names(i) .eq. "NC12H26") iNC12H26 = i
    enddo
    if (iH2.le.0 .or. iH2.gt.nspecies .or. &
         iO2.le.0 .or. iO2.gt.nspecies .or. &
         iCH4.le.0 .or. iCH4.gt.nspecies .or. &
         iN2.le.0 .or. iN2.gt.nspecies .or. &
         iNC12H26.le.0 .or. iNC12H26.gt.nspecies) then
       call bl_pd_abort('Could not find all special species in mech')
    endif

    call bl_pd_is_ioproc(isioproc)

    if (init.ne.1) then
       call bl_pd_abort('probinit called with init ne 1')
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
        
    ! Load domain dimensions into common
    do i=1,amrex_spacedim
       domnlo(i) = problo(i)
       domnhi(i) = probhi(i)
    enddo

    ! Set defaults
    vorterr = 1.e20
    temperr = zero
    adverr = 1.e20
    tempgrad  = 50.0d0
    flametracval = 0.0001d0
    max_temp_lev = 0
    max_vort_lev = 0
    max_trac_lev = 100
    traceSpecVal = 1.d-10
    if (max_vort_lev.lt.0) max_vort_lev=max_temp_lev

    pamb = pphys_getP1atm_MKS()

    splitx = 0.5d0 * (domnhi(1) + domnlo(1))
    xfrontw = 0.05d0 * (domnhi(1) - domnlo(1))
    splity = 0.5d0 * (domnhi(2) + domnlo(2))
    yfrontw = 0.05d0 * (domnhi(2) - domnlo(2))
    blobx = 0.5d0 * (domnhi(1) + domnlo(1))
    bloby = 0.5d0 * (domnhi(2) + domnlo(2))
    blobT = T_in
    Tfrontw = xfrontw
    stTh = -1.d0
    fuel_N2_vol_percent = 0.d0

#if defined(BL_DO_FLCT)
    ! add_turb = .FALSE.
    forceInflow = .FALSE.
    numInflPlanesStore = -1
    forceLo = .TRUE.
    forceHi = .FALSE.
    strmwse_dir = amrex_spacedim
    flct_file = ''
    turb_scale = 1.d0
    nCompInFlow = amrex_spacedim
#endif

    ! Note: for setup with no coflow, set Ro=Rf+wallth
    standoff = zero
    pertmag = 0.d0

    if (isioproc .eq. 1) write(6,*)"reading fortin"
    read(untin,fortin)
    if (isioproc .eq. 1) write(6,*)"done reading fortin"
    read(untin,heattransin)

#if defined(BL_DO_FLCT)
    if (isioproc .eq. 1) write(6,*)"reading flctin"
    read(untin,flctin)
    if (isioproc .eq. 1) write(6,*)"done reading flctin"
#endif

#if defined(BL_DO_FLCT)
    if (forceInflow .eqv. .FALSE.) then
       forceLo = .FALSE.
       forceHi = .FALSE.
    endif
#endif

    ! Set up boundary functions
    if (isioproc .eq. 1) write(6,*)" setup bc"
    call setupbc()
    if (isioproc .eq. 1) write(6,*)" done setup bc"

#if defined(BL_DO_FLCT)
    if (do_flct.eq.1) then
       call init_turbinflow(flct_in, .false.)
    endif
#endif

    if (isioproc.eq.1) then
       write(6,fortin)
       write(6,heattransin)
#if defined(BL_DO_FLCT)
       write(6,flctin)
#endif
    end if

    close(unit=untin)

  end subroutine amrex_probinit
  
  !------------------------------------

  subroutine set_Y_from_Phi(phi,Yt)
    
    use network,   only: nspecies
    use mod_Fvar_def, only : maxspec
    use probdata_module, only : iN2, iO2, iCH4, iH2, H2_frac
    use fuego_chemistry, only : CKXTY

    implicit none

    REAL_T phi, Yt(*)
    
    REAL_T a, alpha,beta,gamma,delt,factor
    REAL_T Xt(maxspec)
    integer n, len

    do n = 1,nspecies
       Xt(n) = 0.d0
    end do
    alpha = H2_frac
    beta = 1.d0 - H2_frac
    gamma = (0.5d0*alpha + 2.d0*beta) / phi
    delt = gamma*.79d0/.21d0
    factor = alpha+beta +gamma+delt
    Xt(iH2) = alpha / factor
    Xt(iCH4) = beta / factor
    Xt(iO2) = gamma / factor
    Xt(iN2) = delt / factor
    CALL CKXTY (Xt, Yt)
  end subroutine set_Y_from_Phi
  
  subroutine setupbc() bind(C, name="setupbc")
    use network,   only: nspecies
    use PeleLM_F,  only: pphys_getP1atm_MKS
    use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : pamb, domnlo, domnhi, maxspec, maxspnml, V_in,&
         fuelID, oxidID, bathID
    use probdata_module, only : Y_bc, T_bc, u_bc, v_bc, w_bc, rho_bc, h_bc, &
         bcinit, T_in, V_co, T_co, fuel_N2_vol_percent, Nzones, iN2, iO2, iH2, iNC12H26
    use user_defined_fcts_3d_module, only : getZone
    use fuego_chemistry, only : CKXTY

    implicit none

    REAL_T Patm, pmf_vals(maxspec+3), a, fact
    REAL_T Xt(maxspec), Yt(maxspec), loc, YF(maxspec), YO(maxspec)
    integer zone, n, fuelZone, airZone, region, len
    integer b(amrex_spacedim),i, j
    data  b / 1, 1, 1 /

    Patm = pamb / pphys_getP1atm_MKS()
    
    ! A diffusion flame
    fuelZone = getZone(domnlo(1), 0.5*(domnlo(2)+domnhi(2)), domnlo(3))
    airZone  = getZone(domnhi(1), domnhi(2), domnhi(3))

    if(fuel_N2_vol_percent .ge. 0.d0)then

       ! Fuel
       do n = 1,nspecies
          Xt(n) = 0.d0
       end do

       Xt(iN2) = fuel_N2_vol_percent*1.d-2
       Xt(iNC12H26) = 1.d0 - Xt(iN2)            

       CALL CKXTY (Xt, Yt)

       do n=1,nspecies
          Y_bc(n-1,fuelZone) = Yt(n)
       end do
       T_bc(fuelZone) = T_in
       u_bc(fuelZone) = 0.d0
       v_bc(fuelZone) = 0.d0
       w_bc(fuelZone) = V_in
          
       ! Air
       do n=1,nspecies
          Xt(n) = zero
       enddo
       Xt(oxidID) = 0.21d0
       Xt(iN2)    = 1.d0 - Xt(oxidID)

       CALL CKXTY (Xt, Yt)
       do n=1,nspecies
          Y_bc(n-1,airZone) = Yt(n)
       end do
          
       T_bc(airZone) = T_co
       u_bc(airZone) = 0.d0
       v_bc(airZone) = 0.d0
       w_bc(airZone) = V_co

    else
       ! Fuel
       do n = 1,nspecies
          Xt(n) = 0.d0
          Yt(n) = 0.d0
       end do

       ! Xt(iNC12H26) = .1208d0
       ! Xt(iO2) = .1319d0
       ! Xt(iN2) = 1.d0 -  Xt(iNC12H26) -  Xt(iO2)

       ! CALL CKXTY (Xt, Yt)

       Yt(iNC12H26) = .3299d0
       Yt(iO2) = .1125d0
       Yt(iN2) = 1.d0 -  Yt(iNC12H26) -  Yt(iO2)

       do n=1,nspecies
          Y_bc(n-1,fuelZone) = Yt(n)
       end do
       T_bc(fuelZone) = T_in
       u_bc(fuelZone) = 0.d0
       v_bc(fuelZone) = 0.d0
       w_bc(fuelZone) = V_in

       ! Air
       do n=1,nspecies
          Xt(n) = zero
       enddo
       Xt(oxidID) = 0.15d0
       Xt(iN2)    = 1.d0 - Xt(oxidID)

       CALL CKXTY (Xt, Yt)         
       do n=1,nspecies
          Y_bc(n-1,airZone) = Yt(n)
       end do
         
       T_bc(airZone) = T_co
       u_bc(airZone) = 0.d0
       v_bc(airZone) = 0.d0
       w_bc(airZone) = V_co
       
    endif

#if 0
    ! do stuff for mixture fraction calculation
    ! read this from input later
    YF=0
    YO=0
    YF(iNC12H26)=1.0
    YO(oxidID) = Y_bc(oxidId-1,airZone)
    YO(iN2) = 1-YO(oxidId)
    Z_fu=0
    Z_ox=0
    do i=1,nspecies
       fact=0
       do j = 1, Nelt
          fact = fact+beta_mix(j)*coeff_mix(i,j)
       enddo
       Z_fu = Z_fu + fact*YF(i)
    enddo
    do i=1,nspecies
       fact=0
       do j = 1, Nelt
          fact = fact+beta_mix(j)*coeff_mix(i,j)
       enddo
       Z_ox = Z_ox + fact*YO(i)
    enddo
#endif

    do zone=1,Nzones
       if (zone.eq.fuelZone .or. zone.eq.airZone) then
         call pphys_RHOfromPTY(b, b, & 
                               rho_bc(zone), DIMARG(b), DIMARG(b), &
                               T_bc(zone),   DIMARG(b), DIMARG(b), &
                               Y_bc(0,zone), DIMARG(b), DIMARG(b), Patm)

         call pphys_HMIXfromTY(b, b, &
                              h_bc(zone),   DIMARG(b), DIMARG(b), &
                              T_bc(zone),   DIMARG(b), DIMARG(b), &
                              Y_bc(0,zone), DIMARG(b), DIMARG(b))

      endif
   enddo

   bcinit = .true.

 end subroutine setupbc

 subroutine init_data(level,time,lo,hi,nscal, &
                      vel,scal,DIMS(state),press,DIMS(press), &
                      delta,xlo,xhi) &
                      bind(C, name="init_data")

   use network,   only: nspecies
   use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
   use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
   use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac
   use mod_Fvar_def, only : domnlo, domnhi, maxspec, maxspnml
   use probdata_module, only : Y_bc, T_co, blobr, Tfrontw, BL_COFLOW
   use user_defined_fcts_3d_module, only : getZone

   implicit none

   integer    level, nscal
   integer    lo(amrex_spacedim), hi(amrex_spacedim)
   integer    DIMDEC(state)
   integer    DIMDEC(press)
   REAL_T     xlo(amrex_spacedim), xhi(amrex_spacedim)
   REAL_T     time, delta(amrex_spacedim)
   REAL_T     vel(DIMV(state),amrex_spacedim)
   REAL_T    scal(DIMV(state),nscal)
   REAL_T   press(DIMV(press))

   integer i, j, k, n, airZone, fuelZone, zone
   integer iO2,iH2,iCH4
   character*(maxspnml) name
   REAL_T x, y, z, ztemp, r, Yl(maxspec), Xl(maxspec), Patm
   REAL_T Xlin(maxspec),alpha,beta,gamma,delt,factor
   REAL_T pmf_vals(maxspec+3), z1, z2, dx, Ly
   REAL_T pert,Lx,FORT_P1ATMMKS,eta,u,v,w,rho,T,h

   fuelZone = getZone(domnlo(1), domnlo(2), domnlo(3))
   airZone  = getZone(domnhi(1), domnhi(2), domnhi(3))

   do k = lo(3), hi(3)
      z = (float(k)+.5)*delta(3)+domnlo(3)
      eta = 0.5d0*(1.d0 - TANH(2.d0*(z-blobr)/Tfrontw))
      do j = lo(2), hi(2)
         y = (float(j)+.5)*delta(2)+domnlo(2)
         do i = lo(1), hi(1)
            x = (float(i)+.5)*delta(1)+domnlo(1)
            do n = 1,nspecies
               scal(i,j,k,FirstSpec+n-1) = Y_bc(n-1,BL_COFLOW)
            end do
            scal(i,j,k,Trac) = 0.d0
            scal(i,j,k,Temp) = T_co
            vel(i,j,k,1) = 0.d0
            vel(i,j,k,2) = 0.d0
            vel(i,j,k,3) = 0.d0
         enddo
      enddo
   enddo

   Patm = pamb / pphys_getP1atm_MKS()

   call pphys_RHOfromPTY(lo,hi, &
        scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state), &
        scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state), &
        scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state), &
        Patm)

   call pphys_HMIXfromTY(lo,hi, &
        scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),RhoH),     DIMS(state), &
        scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state), &
        scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state))

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
 
end module prob_3D_module



