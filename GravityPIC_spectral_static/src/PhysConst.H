#ifndef PHYSCONST_H_
#define PHYSCONST_H_

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <numbers>

namespace PhysConst
{
    inline amrex::Real unit_length_in_cm = 0.0;
    inline amrex::Real unit_vel_in_cm_s  = 0.0;
    inline amrex::Real unit_mass_in_g  = 0.0;
    inline amrex::Real unit_time_in_s = 0.0;

    static constexpr amrex::Real G_cgs = 6.672e-8;
    static constexpr amrex::Real Mpc_in_cm = 3.086e+24;
    static constexpr amrex::Real km_s_in_cm_s = 1e5;

    inline amrex::Real G_code = 0.0;
    inline amrex::Real FourPiG = 0.0;

    inline void Initialize()
    {
        amrex::ParmParse pp;   
        pp.get("unit_length_cm", unit_length_in_cm);
        pp.get("unit_mass_g", unit_mass_in_g);
        pp.get("unit_vel_cm_s",  unit_vel_in_cm_s);
        // convert from cgs to chosen unit system
        unit_time_in_s = unit_length_in_cm / unit_vel_in_cm_s;
        G_code = G_cgs / pow(unit_length_in_cm, 3) * unit_mass_in_g * pow(unit_time_in_s, 2);
        amrex::Print() << "G code : " << G_code << std::endl;
        FourPiG = 4 * std::numbers::pi * G_code;
    }

}
#endif