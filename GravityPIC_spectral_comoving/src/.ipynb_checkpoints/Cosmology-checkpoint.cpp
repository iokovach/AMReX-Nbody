// Cosmology.cpp
#include "Cosmology.H"
#include "PhysConst.H"
#include <AMReX_ParmParse.H>
#include <boost/math/quadrature/gauss_kronrod.hpp>

Cosmology::Cosmology ()
{
    amrex::ParmParse pp;
    pp.get("Omega_m", Omega_m);
    pp.get("Omega_r", Omega_rad);
}

amrex::Real HubbleParam (amrex::Real a, const Cosmology& cosmo) noexcept
{
    return PhysConst::H_0 * std::sqrt(cosmo.Omega_m/(a*a*a)
                               + cosmo.Omega_rad/(a*a*a*a));
}

amrex::Real a_of_t(amrex::Real t) noexcept
{
    amrex::Real t_0 = 2.0 / (3.0*PhysConst::H_0);
    return std::pow(t/t_0, 2.0/3.0);
}

amrex::Real compute_kick_factor (amrex::Real a_old, amrex::Real da, const Cosmology& cosmo)
{
    //using boost::math::quadrature::gauss_kronrod;
    //auto integrand = [&](amrex::Real a){ return 1.0/( a * HubbleParam(a, cosmo)); };
    //return gauss_kronrod<amrex::Real, 15>::integrate(integrand, a_old, a_old+d_a);
    return  2 * (std::sqrt(a_old+da) - std::sqrt(a_old)) / PhysConst::H_0 ;
    //return 2 * (std::sqrt(a_of_t(t_old+dt)) - std::sqrt(a_of_t(t_old))) / PhysConst::H_0 ;
}

amrex::Real compute_drift_factor (amrex::Real a_old, amrex::Real da,
                                          const Cosmology& cosmo)
{
    //using boost::math::quadrature::gauss_kronrod;
    //auto integrand = [&](amrex::Real a){ return 1.0/( a * a * a * HubbleParam(a, cosmo)); };
    //return gauss_kronrod<amrex::Real, 15>::integrate(integrand, a_old, a_old+d_a);
    return  2 * ((1 / std::sqrt(a_old)) - (1 / std::sqrt(a_old+da))) / PhysConst::H_0 ;
    //return  2 * ((1 / std::sqrt(a_of_t(t_old))) - (1 / std::sqrt(a_of_t(t_old+dt)))) / PhysConst::H_0 ;
}