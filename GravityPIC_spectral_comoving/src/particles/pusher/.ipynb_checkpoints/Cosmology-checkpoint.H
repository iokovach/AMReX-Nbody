// Cosmology.H
#ifndef COSMOLOGY_H_
#define COSMOLOGY_H_
#include <AMReX_REAL.H>

struct Cosmology {
    amrex::Real Omega_m;
    amrex::Real Omega_rad;

    Cosmology (); // reads from ParmParse
};

amrex::Real HubbleParam (amrex::Real a, const Cosmology& cosmo) noexcept;
amrex::Real compute_kick_factor  (amrex::Real ln_a_old, amrex::Real d_ln_a, const Cosmology& cosmo);
amrex::Real compute_drift_factor (amrex::Real ln_a_old, amrex::Real d_ln_a, const Cosmology& cosmo);

#endif