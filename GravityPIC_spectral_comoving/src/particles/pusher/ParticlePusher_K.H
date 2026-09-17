#ifndef PARTICLE_PUSHER_K_H_
#define PARTICLE_PUSHER_K_H_

#include <AMReX_REAL.H>

/**
 *
 * \brief This routine advances the particle positions using the given
 * velocity. This is used by the "push_leapfrog" routine below, and also directly to
 * desynchronize the particle positions from the velocities after particle initialization.
 *
 * Arguments:
 *     x         : the particle x position
 *     y         : the particle y position
 *     z         : the particle z position
 *     vx        : the particle x-velocity
 *     vy        : the particle y-velocity
 *     vz        : the particle z-velocity
 *     dt        : the time step
 *     plo       : the LH corner of the problem domain in physical coordinates
 *     phi       : the RH corner of the problem domain in physical coordinates
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void push_leapfrog_positions (amrex::ParticleReal& x,  amrex::ParticleReal& y,  amrex::ParticleReal& z,
                              amrex::ParticleReal& vx, amrex::ParticleReal& vy, amrex::ParticleReal& vz,
                              const amrex::Real drift_factor,
                              const amrex::GpuArray<amrex::Real, 3>& plo,
                              const amrex::GpuArray<amrex::Real, 3>& phi) noexcept
{
    // x_new = x_old + dt * v.
    // (drift_factor = int da/(a^3 H(a)))
    // v here is comoving momentum sans mass a^2 xdot
    x += drift_factor * vx;
    y += drift_factor * vy;
    z += drift_factor * vz;

}

/**
*
* \brief This routine updates the particle positions and velocities using the
* leapfrog time integration algorithm, given the field at the
* particle positions. It also enforces specular reflection off the domain
* walls.
*
* Arguments:
*     x, y, z     : the particle position
*     vx, vy, vz  : the particle velocity
*     Ex, Ey, Ez  : the field at the particle position
*     charge      : the charge of this particle species
*     mass        : the mass of this particle species
*     dt          : the time step
*     plo         : the LH corner of the problem domain in physical coordinates
*     phi         : the RH corner of the problem domain in physical coordinates
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void push_leapfrog (amrex::ParticleReal& x,  amrex::ParticleReal& y,  amrex::ParticleReal& z,
                    amrex::ParticleReal& vx, amrex::ParticleReal& vy, amrex::ParticleReal& vz,
                    const amrex::ParticleReal Ex, const amrex::ParticleReal Ey, const amrex::ParticleReal Ez,
                    const amrex::Real kick_factor, const amrex::Real drift_factor, const amrex::GpuArray<amrex::Real, 3>& plo,
                    const amrex::GpuArray<amrex::Real, 3>& phi) noexcept
{
    // (kick_factor = int da/ a^2 H(a))
    // Ex = - grad phi 
    vx += kick_factor * Ex;
    vy += kick_factor * Ey;
    vz += kick_factor * Ez;

    // update positions using the newly-updated velocities.
    push_leapfrog_positions(x, y, z, vx, vy, vz, drift_factor, plo, phi);
}

#endif