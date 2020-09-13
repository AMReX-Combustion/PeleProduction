#include <extern_parameters.H>
#include <extern_parameters_F.H>

  AMREX_GPU_MANAGED amrex::Real small_massfrac;

  AMREX_GPU_MANAGED amrex::Real mwt_scalar;

  AMREX_GPU_MANAGED int numaux;

  std::string auxnamesin;

  AMREX_GPU_MANAGED amrex::Real react_T_min;

  AMREX_GPU_MANAGED amrex::Real react_T_max;

  AMREX_GPU_MANAGED amrex::Real react_rho_min;

  AMREX_GPU_MANAGED amrex::Real react_rho_max;

  AMREX_GPU_MANAGED int new_Jacobian_each_cell;

  AMREX_GPU_MANAGED int use_bulk_viscosity;

  AMREX_GPU_MANAGED int do_flct;

  AMREX_GPU_MANAGED int flct_dir;

  std::string flct_in;


  void init_extern_parameters() {
    int slen = 0;

    get_f90_small_massfrac(&small_massfrac);

    get_f90_mwt_scalar(&mwt_scalar);

    get_f90_numaux(&numaux);

    get_f90_auxnamesin_len(slen);
    char _auxnamesin[slen+1];
    get_f90_auxnamesin(_auxnamesin);
    auxnamesin = std::string(_auxnamesin);

    get_f90_react_T_min(&react_T_min);

    get_f90_react_T_max(&react_T_max);

    get_f90_react_rho_min(&react_rho_min);

    get_f90_react_rho_max(&react_rho_max);

    get_f90_new_Jacobian_each_cell(&new_Jacobian_each_cell);

    get_f90_use_bulk_viscosity(&use_bulk_viscosity);

    get_f90_do_flct(&do_flct);

    get_f90_flct_dir(&flct_dir);

    get_f90_flct_in_len(slen);
    char _flct_in[slen+1];
    get_f90_flct_in(_flct_in);
    flct_in = std::string(_flct_in);

  }
