#include <PeleLM.H>
#include <pelelm_prob.H>
#include <pmf.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        amrex::ParmParse pp("prob");

        pp.query("P_mean",   PeleLM::prob_parm->P_mean);
        pp.query("standoff", PeleLM::prob_parm->standoff);
        pp.query("pertmag",  PeleLM::prob_parm->pertmag);
        pp.query("Vin", PeleLM::prob_parm->Vin);
        pp.query("Vcoflow", PeleLM::prob_parm->Vcoflow);
        pp.query("slot_width", PeleLM::prob_parm->slot_width);
        pp.query("Tstart", PeleLM::prob_parm->Tstart);

    }
}
