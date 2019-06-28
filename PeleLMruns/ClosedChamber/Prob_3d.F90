#include <AMReX_REAL.H>
#include <AMReX_CONSTANTS.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ArrayLim.H>

#include <Prob_F.H>
#include <PeleLM_F.H>


module prob_3D_module

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
  
  subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(C, name="amrex_probinit")

    use PeleLM_F,  only: pphys_getP1atm_MKS
    use mod_Fvar_def, only : pamb, dpdt_factor, closed_chamber
    use mod_Fvar_def, only : fuelID, domnhi, domnlo, dim
    use probdata_module, only : standoff, pertmag, pert_scale, rho_bc, Y_bc, blobz

    implicit none

    integer init, namlen
    integer name(namlen)
    integer untin
    REAL_T problo(dim), probhi(dim)

    integer i
    REAL_T area

    namelist /fortin/ standoff, pertmag, pert_scale, blobz
    namelist /heattransin/ pamb, dpdt_factor, closed_chamber

    !
    !      Build `probin' filename -- the name of file containing fortin namelist.
    !
    integer maxlen, isioproc
    parameter (maxlen=256)
    character probin*(maxlen)

    call bl_pd_is_ioproc(isioproc)

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

    !     Load domain dimensions into common (put something computable there for SDIM<3)
    domnlo = problo
    domnhi = probhi

    untin = 9
    open(untin,file=probin(1:namlen),form='formatted',status='old')

    pamb = pphys_getP1atm_MKS()
    dpdt_factor = 1.0
    closed_chamber = 0

    !     Note: for setup with no coflow, set Ro=Rf+wallth
    standoff = zero
    pertmag = 0.d0
    pert_scale = 1.0d0

    if (isioproc .eq. 1) then
       write(6,*)"reading fortin"
    endif

    read(untin,fortin)

    if (isioproc .eq. 1) then
       write(6,*)"done reading fortin"
    endif

    read(untin,heattransin)

    close(unit=untin)

    !     Set up boundary functions
    call setupbc()

    if (isioproc.eq.1) then
       write(6,fortin)
       write(6,heattransin)
    end if

  end subroutine amrex_probinit

  subroutine setupbc() bind(C, name="setupbc")
    
    use network,   only: nspec
    use PeleLM_F,  only: pphys_getP1atm_MKS
    use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use fuego_chemistry, only: ckxty
    use mod_Fvar_def, only : pamb
    use probdata_module, only : standoff, rho_bc, u_bc, v_bc, w_bc, T_bc, h_bc, Y_bc, bcinit, domnlo

    implicit none

    REAL_T Patm, pmf_vals(nspec+3), a
    REAL_T Xt(nspec), loc
    integer n, len, b(3)

    b = 1
    Patm = pamb / pphys_getP1atm_MKS()

    ! Take fuel mixture from pmf file
    loc = (domnlo(2)-standoff)*100.d0
    call pmf(loc,loc,pmf_vals,n)

    if (n.ne.nspec+3) then
       call bl_pd_abort('INITDATA: n(pmf) .ne. Nspec+3')
    endif

    do n = 1,Nspec
       Xt(n) = pmf_vals(3+n)
    end do

    CALL ckxty (Xt, Y_bc(0))

    T_bc(1) = pmf_vals(1)
    u_bc = zero
    v_bc = zero
    w_bc = zero
    
    ! Set density and hmix consistent with data
    call pphys_RHOfromPTY(b, b, &
                          rho_bc(1), DIMARG(b), DIMARG(b), &
                          T_bc(1),   DIMARG(b), DIMARG(b), &
                          Y_bc(0),   DIMARG(b), DIMARG(b), Patm)
    call pphys_HMIXfromTY(b, b, &
                          h_bc(1),   DIMARG(b), DIMARG(b), &
                          T_bc(1),   DIMARG(b), DIMARG(b), &
                          Y_bc(0),   DIMARG(b), DIMARG(b))

    bcinit = .true.

  end subroutine setupbc

  ! ::: -----------------------------------------------------------

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
       delta,xlo,xhi) bind(C, name="init_data")
    
    use network,   only: nspec
    use PeleLM_F,  only: pphys_getP1atm_MKS, pphys_get_spec_name2
    use PeleLM_3D, only: pphys_RHOfromPTY, pphys_HMIXfromTY
    use fuego_chemistry, only: ckxty
    use mod_Fvar_def, only : Density, Temp, FirstSpec, RhoH, pamb, Trac, dim
    use probdata_module, only : standoff, pert_scale, pertmag, domnlo, domnhi, blobz

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

    integer i, j, k, n
    REAL_T x, y, z, r, Yl(nspec), Xl(nspec), Patm
    REAL_T pmf_vals(nspec+3), z1, z2, dx, Ly
    REAL_T pert,Lx,eta,u,v,w,rho,T,h

    do k = lo(3), hi(3)
       z = (float(k)+.5d0)*delta(3)+domnlo(3)
       do j = lo(2), hi(2)
          y = (float(j)+.5d0)*delta(2)+domnlo(2)
          do i = lo(1), hi(1)
             x = (float(i)+.5d0)*delta(1)+domnlo(1)

             pert = 0.d0
             if (pertmag .gt. 0.d0) then
                Lx = (domnhi(1) - domnlo(1)) / pert_scale
                Ly = (domnhi(2) - domnlo(2)) / pert_scale
                pert = pertmag*(1.000 * sin(2*Pi*4*x/Lx)             * sin(2*Pi*5*y/Ly)&
                     + 1.023 * sin(2*Pi*2*(x-.004598)/Lx)   * sin(2*Pi*4*(y-.0053765)/Ly)&
                     + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx) * sin(2*Pi*3*(y-.02137)/Ly)&
                     + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)     * sin(2*Pi*6*(y-.018)/Ly)&
                     + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
                
             endif

             z1 = (z - blobz - standoff - 0.5d0*delta(2) + pert )*100.d0
             z2 = (z - blobz - standoff + 0.5d0*delta(2) + pert )*100.d0 

             call pmf(z1,z2,pmf_vals,nPMF)               
             if (nPMF.ne.Nspec+3) then
                call bl_abort('INITDATA: n .ne. Nspec+3')
             endif

             scal(i,j,k,Temp) = pmf_vals(1)
             do n = 1,Nspec
                Xl(n) = pmf_vals(3+n)
             end do

             CALL ckxty (Xl, Yl)

             do n = 1,Nspec
                scal(i,j,k,FirstSpec+n-1) = Yl(n)
             end do

             scal(i,j,k,Trac) = 0.d0

             vel(i,j,k,1) = 0.d0
             vel(i,j,k,2) = 0.d0
             vel(i,j,k,3) = pmf_vals(2)*1.d-2

          end do
       end do
    end do

    Patm = pamb / pphys_getP1atm_MKS()

    call pphys_RHOfromPTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Density),  DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state),&
         Patm)

    call pphys_HMIXfromTY(lo,hi,&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),RhoH),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),Temp),     DIMS(state),&
         scal(ARG_L1(state),ARG_L2(state),ARG_L3(state),FirstSpec),DIMS(state))

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             do n = 0,Nspec-1
                scal(i,j,k,FirstSpec+n) = scal(i,j,k,FirstSpec+n)*scal(i,j,k,Density)
             enddo
             scal(i,j,k,RhoH) = scal(i,j,k,RhoH)*scal(i,j,k,Density)
          enddo
       enddo
    enddo
  end subroutine init_data

end module prob_3D_module
