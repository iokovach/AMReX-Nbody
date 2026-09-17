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
amrex::Real a_of_t (amrex::Real t) noexcept;
amrex::Real compute_kick_factor  (amrex::Real a_old, amrex::Real da, const Cosmology& cosmo);
amrex::Real compute_drift_factor (amrex::Real a_old, amrex::Real da, const Cosmology& cosmo);

#endif