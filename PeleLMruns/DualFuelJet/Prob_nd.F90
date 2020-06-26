
#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>
#include "mechanism.h"

module prob_nD_module

  use amrex_fort_module, only : dim=>amrex_spacedim
  use amrex_error_module, only : amrex_abort
  use amrex_omp_module

  implicit none

  private

  public :: amrex_probinit, setupbc, init_data

  integer, parameter :: nspecies = NUM_SPECIES

contains

  ! ::: -----------------------------------------------------------

  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

#if defined(BL_DO_FLCT)
    use turbinflow_module, only : init_turbinflow
#endif
    use PeleLM_F,  only : pphys_getP1atm_MKS
    !use derive_PLM_3D, only: init_mixture_fraction
    use mod_Fvar_def, only : pamb
    use probdata_module, only : V_in, Y_ox, Y_fu
    use probdata_module, only : V_in, Y_fu, Y_ox, V_co, &
                                T_ox, T_fu, mixfrac_in, mixfrac_co, h_fu_nist, &
                                splitx, xfrontw, &
                                blobr, bloby, blobx, blobz, blobT, Tfrontw, turb_scale, &
                                domnlo, domnhi, adverr, flametracval, H2_frac, &
                                max_diff_lev, max_hrr_lev, max_mix_lev, max_temp_lev, &
                                max_trac_lev, max_vort_lev, pertmag, &
                                pseudo_gravity, refine_nozzle, splity, standoff, &
                                stTh, T_switch, temperr, tempgrad, traceSpecVal, &
                                vorterr, yfrontw, &
#if defined(BL_DO_FLCT)
                                do_flct, &
#endif
                                iN2
    use probdata_module, only: jet_mixfrac, jet_hrr, jet_zmax, &
                               jet_core_mixfrac, jet_core_zmax, &
                               ltc_RO2, &
                               edge_flame_OH, edge_flame_zmax, &
                               prmx_flame_HO2, prmx_flame_mixfrac, prmx_flame_zmax
    use fuego_chemistry, only : get_species_index, spec_names

    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(dim), probhi(dim)

    integer i
    character flct_file*(72)
    character streamsin_file*(72)

    namelist /fortin/ vorterr, temperr, adverr, tempgrad, &
                      flametracval, max_temp_lev, max_vort_lev, &
                      max_trac_lev, max_hrr_lev, max_mix_lev, &
                      traceSpecVal, &
                      turb_scale, V_in, V_co, H2_frac,T_switch, &
                      standoff, pertmag, splitx, xfrontw, &
                      splity, yfrontw, blobx, bloby, blobz, blobr, &
                      blobT, Tfrontw, stTh, refine_nozzle, &
                      max_diff_lev, streamsin_file, &
                      T_ox, T_fu, mixfrac_in, mixfrac_co, h_fu_nist

    namelist /heattransin/ pamb
#if defined(BL_DO_FLCT)
    namelist /flctin/ do_flct, flct_file, turb_scale
#endif

    namelist /CMA/ jet_mixfrac, jet_hrr, jet_zmax, &
                   jet_core_mixfrac, jet_core_zmax, &
                   ltc_RO2, edge_flame_OH, edge_flame_zmax, &
                   prmx_flame_HO2, prmx_flame_mixfrac, prmx_flame_zmax

    !
    ! Build `probin' filename -- the name of file containing fortin namelist.
    !
    integer maxlen, isioproc
    parameter (maxlen=256)
    character probin*(maxlen)

    do i = 1,nspecies
       if (spec_names(i) .eq. "N2") iN2 = i
    enddo
    if (iN2.lt.1 .or. iN2.gt.nspecies) then
       call bl_pd_abort('Could not find N2 in the mechanism')
    endif

    call bl_pd_is_ioproc(isioproc)

!    if (init.ne.1) then
!       call bl_pd_abort('probinit called with init ne 1')
!    end if

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
    do i=1,dim
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
    blobT = T_fu
    Tfrontw = xfrontw
    stTh = -1.d0

    ! IC/BC
    mixfrac_in = -1.0
    mixfrac_co = zero
    T_fu = -1.0
    T_ox = -1.0
    h_fu_nist = .true.
    streamsin_file = 'streams.in'

#if defined(BL_DO_FLCT)
    do_flct = .FALSE.
    flct_file = ''
#endif

    ! Note: for setup with no coflow, set Ro=Rf+wallth
    standoff = zero
    pertmag = zero

    if (isioproc .eq. 1) write(6,*)"reading fortin"
    read(untin,fortin)
    if (isioproc .eq. 1) write(6,*)"done reading fortin"

    if (isioproc .eq. 1) write(6,*)"reading heattransin"
    read(untin,heattransin)
    if (isioproc .eq. 1) write(6,*)"done reading heattransin"

#if defined(BL_DO_FLCT)
    if (isioproc .eq. 1) write(6,*)"reading flctin"
    read(untin,flctin)
    if (isioproc .eq. 1) write(6,*)"done reading flctin"
#endif
    if (isioproc .eq. 1) write(6,*)"reading CMA"
    read(untin,CMA)
    if (isioproc .eq. 1) write(6,*)"done reading CMA"

    ! Close probin
    close(unit=untin)

    ! Write probin data to screen
    if (isioproc.eq.1) then
       write(6,fortin)
       write(6,heattransin)
#if defined(BL_DO_FLCT)
       write(6,flctin)
#endif
       write(6,CMA)
    end if

    if (isioproc .eq. 1) write(6,*)" setup bc"
    ! Initialise the mixture fraction data
    call read_streams_input(streamsin_file)
    !call init_mixture_fraction()

    ! Set up boundary functions
    call setupbc()
    if (isioproc .eq. 1) write(6,*)" done setup bc"

#if defined(BL_DO_FLCT)
    ! Initialise turbulent inflow data
    if (do_flct) then
       !$omp parallel
       call init_turbinflow(flct_file, .false.)
       !$omp end parallel
    endif
#endif

  end subroutine amrex_probinit

  ! ::: -----------------------------------------------------------

  subroutine read_streams_input(streamsin)

    use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor, &
                                                amrex_pd_bcast
    use amrex_string_module, only : amrex_string_f_to_lower
    use probdata_module, only : Y_fu, Y_ox
    use fuego_chemistry
    
    implicit none

! In/Out
    character*(*), intent(in) :: streamsin

! Local
    REAL_T, dimension(nspecies) :: X_ox, X_fu
    logical, dimension(nspecies) :: species_nzero
    logical :: is_ioproc, exists
    integer, parameter :: infile = 7
    REAL_T :: tmp
    integer :: stat, i, is_mole
    character(len=1024) :: buf
    character(len=64) :: name, value
    REAL_T, parameter :: tol = ten*tiny(zero)


    is_ioproc = amrex_pd_ioprocessor()
    is_mole = -1

    ! Check if streams input file exists.
    ! NOTE: AMReX doesn't have a Fortran broadcast function for logical or
    !       integer variables
    tmp = -1.0
    if (is_ioproc) then
      inquire(file=streamsin, exist=exists)
      if (exists) tmp = 1.0
    endif
    call amrex_pd_bcast(tmp)
    if (tmp < 0.0) call amrex_abort('streams input file does not exist!')

    ! Only the I/O processor reads the file
    if (is_ioproc) then
      ! Set mole fractions initially to zero
      X_ox = 0.0
      X_fu = 0.0

      ! Read and parse the file
      open(unit=infile, file=streamsin, form='formatted', status='old' )
      do
        read(infile, '(A)', iostat=stat) buf
        if (stat .ne. 0) exit

        ! Read keywords to determine what to do (allow this to be case insensitive)
        buf = adjustl(amrex_string_f_to_lower(buf))
        if ((index(buf, '!') .eq. 1) .or. (index(buf, '/') .eq. 1)) then
          ! Skip comment lines
          cycle
        else if (index(buf, 'mole') .eq. 1) then
          if (is_mole < 0) then
            is_mole = 1
          else
            call amrex_abort('Multiple "mass" or "mole" keywords!')
          endif
        else if (index(buf, 'mass') .eq. 1) then
          if (is_mole < 0) then
            is_mole = 1
          else
            call amrex_abort('Multiple "mass" or "mole" keywords!')
          endif
        else if (index(buf, 'fuel') .eq. 1) then
          do
            read(infile, '(A)', iostat=stat) buf
            if (index(amrex_string_f_to_lower(buf), 'end') .gt. 0) exit
            call get_name_value(name, value, buf)
            i = get_species_index(name)
            read(value,*) X_fu(i)
          end do
        else if (index(buf, 'oxid') .eq. 1) then
          do
            read(infile, '(A)', iostat=stat) buf
            if (index(amrex_string_f_to_lower(buf), 'end') .gt. 0) exit
            call get_name_value(name, value, buf)
            i = get_species_index(name)
            read(value,*) X_ox(i)
          end do
        else
          ! Empty lines
          cycle
        endif
      end do
      close(unit=infile)

      ! Convert to mass fractions if necessary
      if (is_mole .eq. 0) then
        Y_fu = X_fu
        Y_ox = X_ox
      else if (is_mole .eq. 1) then
        call ckxty(X_fu,Y_fu)
        call ckxty(X_ox,Y_ox)
      else
        call amrex_abort('Could not determine if mole or mass fractions!')
      endif
    endif

    ! Broadcast the data to all other processors
    call amrex_pd_bcast(Y_fu)
    call amrex_pd_bcast(Y_ox)

    ! Normalise compositions if necessary
    tmp = sum(Y_fu)
    if ((tmp .gt. one) .or. (tmp .lt. one-tol)) then
      if (is_ioproc) write(6,*) 'WARNING: Normalizing fuel composition...'
      Y_fu = Y_fu/tmp
    endif

    tmp = sum(Y_ox)
    if ((tmp .gt. one) .or. (tmp .lt. one-tol)) then
      if (is_ioproc) write(6,*) 'WARNING: Normalizing oxidiser composition...'
      Y_ox = Y_ox/tmp
    endif

    contains
      subroutine get_name_value(name, value, str)
        implicit none

        character*(*), intent(inout) :: str
        character*(*), intent(out) :: name, value
        integer :: i

        str = adjustl(str)

        i = index(str, ' ')

        name = adjustl(str(:i-1))
        value = adjustl(str(i+1:))
      end subroutine get_name_value
  end subroutine read_streams_input

  ! ::: -----------------------------------------------------------

  subroutine setupbc() bind(C, name="setupbc")
    use PeleLM_nd, only: pphys_RHOfromPTY, pphys_HMIXfromTY, pphys_TfromHY
    use PeleLM_F, only: pphys_getP1atm_MKS
    use amrex_paralleldescriptor_module, only : amrex_pd_ioprocessor
    use mod_Fvar_def, only : pamb
    use probdata_module, only : V_in, Y_ox, Y_fu
    use probdata_module, only : bcinit, BL_FUELPIPE, BL_COFLOW, &
                                Y_bc, T_bc, u_bc, v_bc, w_bc, rho_bc, h_bc, &
                                V_co, Nzones, iN2, h_fu_nist, &
                                T_ox, T_fu, mixfrac_in, mixfrac_co
    use probdata_module, only : get_nist_enthalpy, get_T_from_hY
    use fuego_chemistry, only : spec_names

    implicit none

    REAL_T Patm, Yt, h_ox, h_fu
    integer i, j
    logical :: is_ioproc
    integer :: n
    integer :: Niter
    REAL_T,parameter :: Ttol = 1.d-12
    integer, parameter :: NiterMAX = 20
    REAL_T :: Tres(NiterMAX)
    integer,parameter,dimension(3) :: b_lo = (/1,1,1/)
    integer,parameter,dimension(3) :: b_hi = (/1,1,1/)
    REAL_T harr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1))
    REAL_T Tarr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1))
    REAL_T Yarr(b_lo(1):b_hi(1),b_lo(1):b_hi(1),b_lo(1):b_hi(1),nspecies)

    is_ioproc = amrex_pd_ioprocessor()


    if (     (mixfrac_in .lt. zero) .and. (mixfrac_in .gt. one) &
        .or. (mixfrac_co .lt. zero) .and. (mixfrac_co .gt. one) &
        .or.                                   (T_fu .lt. zero) &
        .or.                                   (T_ox .lt. zero)) then
      call amrex_abort('Invalid mixfrac_in, mixfrac_co, T_ox or T_fu')
    endif


    ! --------------------------------------------------------------------------
    ! Enthalpy in the pure fuel and oxidiser streams
    ! --------------------------------------------------------------------------
    ! Oxidiser enthalpy
    Tarr = T_ox
    Yarr(1,1,1,:) = Y_ox(:)
    call pphys_HMIXfromTY(b_lo, b_hi, &
                          harr, b_lo, b_hi, &
                          Tarr, b_lo, b_hi, &
                          Yarr, b_lo, b_hi)
    h_ox = harr(1,1,1)

    ! Fuel enthalpy - here we account for the heat of evaporation
    Tarr = T_fu
    Yarr(1,1,1,:) = Y_fu(:)
    call pphys_HMIXfromTY(b_lo, b_hi, &
                          harr, b_lo, b_hi, &
                          Tarr, b_lo, b_hi, &
                          Yarr, b_lo, b_hi)
    h_fu = harr(1,1,1)
    
    if (h_fu_nist) then
      if(is_ioproc) then
        write(6,*)
        write(6,'(2x,a)') 'Apply enthalpy correction'
        write(6,'(4x,a,1x,es12.5,1x,a)') 'Fuel enthalpy prior correction: h =', &
                                          h_fu, 'J/kg'
      endif

      h_fu = 0.0
      do i=1,nspecies
        ! Only consider species that are present in the fuel stream
        if(Y_fu(i) .gt. 1d-12) then
          h_fu = h_fu + Y_fu(i)*get_nist_enthalpy(spec_names(i), T_fu)
        endif
      enddo

      if(is_ioproc) then
        write(6,'(4x,a,1x,es12.5,1x,a)') 'Fuel enthalpy after correction: h =', &
                                          h_fu, 'J/kg'
      endif
    endif


    ! --------------------------------------------------------------------------
    ! Jet (fuel-pipe)
    ! --------------------------------------------------------------------------
    ! Mass fractions and enthalpy are a linear combination of the fuel and
    ! oxidiser streams
    Yt = zero
    do i=1,nspecies
      if (i .ne. iN2) then
        Y_bc(i-1,BL_FUELPIPE) = Y_ox(i) + (Y_fu(i) - Y_ox(i))*mixfrac_in
        Yt = Yt + Y_bc(i-1,BL_FUELPIPE)
      endif
    enddo
    Y_bc(iN2-1,BL_FUELPIPE) = one - Yt
    h_bc(BL_FUELPIPE) = h_ox + (h_fu - h_ox)*mixfrac_in

    ! Get temperature from Y and h; linearly interpolate the initial guess
    T_bc(BL_FUELPIPE) = T_ox + (T_fu - T_ox)*mixfrac_in
    Niter = pphys_TfromHY(b_lo, b_hi, &
                          T_bc(BL_FUELPIPE),              b_lo, b_hi, &
                          h_bc(BL_FUELPIPE),              b_lo, b_hi, &
                          Y_bc(0:nspecies-1,BL_FUELPIPE), b_lo, b_hi, &
                          Ttol, NiterMAX, Tres)

    ! Density
    Patm = pamb / pphys_getP1atm_MKS()
    call pphys_RHOfromPTY(b_lo, b_hi, &
                          rho_bc(BL_FUELPIPE),            b_lo, b_hi, &
                          T_bc(BL_FUELPIPE),              b_lo, b_hi, &
                          Y_bc(0:nspecies-1,BL_FUELPIPE), b_lo, b_hi, Patm)

    ! Velocity
    v_bc(BL_FUELPIPE) = zero
    u_bc(BL_FUELPIPE) = zero
    w_bc(BL_FUELPIPE) = V_in


    ! --------------------------------------------------------------------------
    ! Co-flow
    ! --------------------------------------------------------------------------
    ! Mass fractions and enthalpy are a linear combination of the fuel and
    ! oxidiser streams
    Yt = zero
    do i=1,nspecies
      if (i .ne. iN2) then
        Y_bc(i-1,BL_COFLOW) = Y_ox(i) + (Y_fu(i) - Y_ox(i))*mixfrac_co
        Yt = Yt + Y_bc(i-1,BL_COFLOW)
      endif
    enddo
    Y_bc(iN2-1,BL_COFLOW) = one - Yt
    h_bc(BL_COFLOW) = h_ox + (h_fu - h_ox)*mixfrac_co

    ! Get temperature from Y and h; linearly interpolate the initial guess
    T_bc(BL_COFLOW) = T_ox + (T_fu - T_ox)*mixfrac_co
    Niter = pphys_TfromHY(b_lo, b_hi, &
                          T_bc(BL_COFLOW),              b_lo, b_hi, &
                          h_bc(BL_COFLOW),              b_lo, b_hi, &
                          Y_bc(0:nspecies-1,BL_COFLOW), b_lo, b_hi, &
                          Ttol, NiterMAX, Tres)

    ! Density
    call pphys_RHOfromPTY(b_lo, b_hi, &
                          rho_bc(BL_COFLOW),            b_lo, b_hi, &
                          T_bc(BL_COFLOW),              b_lo, b_hi, &
                          Y_bc(0:nspecies-1,BL_COFLOW), b_lo, b_hi, Patm)

    ! Velocity
    u_bc(BL_COFLOW) = zero
    v_bc(BL_COFLOW) = zero
    w_bc(BL_COFLOW) = V_co


    ! --------------------------------------------------------------------------
    ! Print BCs
    ! --------------------------------------------------------------------------
    if (is_ioproc) then
      write(6,'(2x,a)') 'Inlet boundary conditions:'
      do j=1,Nzones
        if (j.eq.BL_FUELPIPE) write(6,'(4x,a)') 'Jet properties:'
        if (j.eq.BL_COFLOW)   write(6,'(4x,a)') 'Co-flow properties:'
        if (j.eq.BL_FUELPIPE .or. j.eq.BL_COFLOW) then
          write(6,'(4x,a22,1x,f12.7)') 'u = ', u_bc(j)
          write(6,'(4x,a22,1x,f12.7)') 'v = ', v_bc(j)
          write(6,'(4x,a22,1x,f12.7)') 'w = ', w_bc(j)
          write(6,'(4x,a22,1x,f12.7)') 'T = ', T_bc(j)
          write(6,'(4x,a22,1x,f12.1)') 'h = ', h_bc(j)
          write(6,'(4x,a22,1x,f12.7)') 'rho = ', rho_bc(j)
          do i = 1, nspecies
            if (Y_bc(i-1,j) .gt. 1e-14) then
              write(6,'(4x,a22,1x,f12.7)') adjustl(spec_names(i)), Y_bc(i-1,j)
            endif
          enddo
        endif
      enddo
    endif

    bcinit = .true.

  end subroutine setupbc

  ! ::: -----------------------------------------------------------

  subroutine init_data(level, time, lo, hi, nscal, &
                       vel, scal, s_lo, s_hi, press, p_lo, p_hi, &
                       delta, xlo, xhi) &
                       bind(C, name="init_data")

    use PeleLM_F,  only: pphys_getP1atm_MKS
    use PeleLM_nd, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb
    use mod_Fvar_def, only : domnlo, domnhi
    use probdata_module, only : Y_bc, T_bc, blobr, Tfrontw, BL_COFLOW, splitx, xfrontw
    use user_defined_fcts_nd_module, only : bcfunction

    implicit none

    integer, intent(in) :: level, nscal
    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: p_lo(3), p_hi(3)
    REAL_T, intent(in)  :: xlo(3), xhi(3)
    REAL_T, intent(in)  :: time, delta(3)
    REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),dim), intent(out) :: vel
    REAL_T, dimension(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),nscal), intent(out) :: scal
    REAL_T, dimension(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3)), intent(out) :: press

    integer i, j, k, n
    integer iO2,iH2,iCH4
    REAL_T x, y, z, ztemp, r, Yl(nspecies), Xl(nspecies), Patm
    REAL_T Xlin(nspecies),alpha,beta,gamma,delt,factor
    REAL_T pmf_vals(nspecies+3), z1, z2, dx, Ly
    REAL_T pert,Lx,FORT_P1ATMMKS,eta,u,v,w,rho,T,h

    REAL_T xv(3),velv(3),TT,hh,eta1

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
          scal(i,j,k,Temp) = T_bc(BL_COFLOW)

          vel(i,j,k,1) = zero
          vel(i,j,k,2) = zero
          vel(i,j,k,3) = zero

          xv(1) = x
          xv(2) = y
          xv(3) = z
          call bcfunction(xv,delta,3,1,time,.true.,velv,rho,Yl,TT,hh)

          eta1 = 0.5d0*(1.d0 - TANH(2.d0*(z-blobr)/Tfrontw))
          do n = 1,nspecies
            scal(i,j,k,FirstSpec+n-1) = eta1*Yl(n) + (1.d0-eta1)*scal(i,j,k,FirstSpec+n-1)
          end do
          scal(i,j,k,Temp) = eta1*TT + (1.d0-eta1)*scal(i,j,k,Temp)

          do n=1,3
            vel(i,j,k,n) = eta1*velv(n) + (1.d0-eta1)*vel(i,j,k,n)
          enddo
          
        enddo
      enddo
    enddo

    Patm = pamb / pphys_getP1atm_MKS()

    call pphys_RHOfromPTY(lo,hi, &
          scal(s_lo(1),s_lo(2),s_lo(3),Density),  s_lo, s_hi, &
          scal(s_lo(1),s_lo(2),s_lo(3),Temp),     s_lo, s_hi, &
          scal(s_lo(1),s_lo(2),s_lo(3),FirstSpec),s_lo, s_hi, &
          Patm)

    call pphys_HMIXfromTY(lo,hi, &
          scal(s_lo(1),s_lo(2),s_lo(3),RhoH),     s_lo, s_hi, &
          scal(s_lo(1),s_lo(2),s_lo(3),Temp),     s_lo, s_hi, &
          scal(s_lo(1),s_lo(2),s_lo(3),FirstSpec),s_lo, s_hi)

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

end module prob_nD_module
