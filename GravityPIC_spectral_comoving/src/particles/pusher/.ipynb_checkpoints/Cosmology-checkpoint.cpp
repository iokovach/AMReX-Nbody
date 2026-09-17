// Cosmology.cpp
#include "Cosmology.H"
#include <AMReX_ParmParse.H>
#include <boost/math/quadrature/gauss_kronrod.hpp>

Cosmology::Cosmology ()
{
    amrex::ParmParse pp;
    pp.get("Omega_m", Omega_m);
    pp.get("Omega_r", Omega_r);
}

amrex::Real HubbleParam (amrex::Real a, const Cosmology& cosmo) noexcept
{
    return PhysConst::H0 * std::sqrt(cosmo.Omega_m/(a*a*a)
                               + cosmo.Omega_r/(a*a*a*a));
}

amrex::Real compute_kick_factor (amrex::Real ln_a_old, amrex::Real d_ln_a, const Cosmology& cosmo)
{
    using boost::math::quadrature::gauss_kronrod;
    auto integrand = [&](amrex::Real a){ return 1.0/(a * a * HubbleParam(a, cosmo)); };
    return gauss_kronrod<amrex::Real, 15>::integrate(integrand, std::exp(ln_a_old), std::exp(ln_a_old+d_ln_a));
}

amrex::Real compute_drift_factor (amrex::Real ln_a_old, amrex::Real d_ln_a,
                                          const Cosmology& cosmo)
{
    using boost::math::quadrature::gauss_kronrod;
    auto integrand = [&](amrex::Real a){ return 1.0/(a * a * a * HubbleParam(a, cosmo)); };
    return gauss_kronrod<amrex::Real, 15>::integrate(integrand, std::exp(ln_a_old), std::exp(ln_a_old+d_ln_a));
}