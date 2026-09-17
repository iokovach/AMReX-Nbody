#include "GravityPICAmr.H"
#include "field_solver/FieldSolver.H"
#include "particles/GravityParticleContainer.H"
#include "diagnostics/FieldIO.H"
#include "Cosmology.H"

#include <cmath>

// constructor -- make an AmrCore object and read in simulation parameters
GravityPICAmr::GravityPICAmr()
    : amrex::AmrCore()
{
    amrex::ParmParse pp;
    amrex::Array<int, AMREX_SPACEDIM> pdc;
    
    pp.query("max_level", max_level);
    pp.query("rho_err", rho_err);
    pp.query("regrid_int", regrid_int);
    pp.get("n_buffer", n_buffer);
    pp.get("max_step", max_step);
    pp.get("dt", dt);
    pp.get("d_ln_a", d_ln_a);
    pp.get("ic_file", ic_file);
    pp.get("geometry.is_periodic", pdc);

    pp.get("a_ini", a_ini);

    pp.query("plot_int", plot_int);
    pp.query("particle_output_int", particle_output_int);

    // get bc as an int for phi solver
    bc = pdc[0];

}

// function to build boxarrays and particles
void
GravityPICAmr::InitData()
{
    // Build the initial BoxArray, Geom and DM at level 0 
    InitFromScratch(0.0);

    // mask out coarse and fine regions for field gathering later
    RebuildMasks();

    // link PC to GravityPICAmr object instead of just box arrays/geom etc at step 0
    myPC = std::make_unique<GravityParticleContainer>(this);
    
    // read particle data and put them in the container
    myPC->InitParticles(ic_file);
    
    amrex::Print() << "InitData: loaded " << myPC->TotalNumberOfParticles()
                   << " particles\n";
}

// detect boundaries between coarse/fine regions
void GravityPICAmr::RebuildMasks()
{
    const int nlev = finestLevel() + 1;

    // Get the current AMR hierarchy from AmrCore
    amrex::Vector<amrex::BoxArray> grids(nlev);
    amrex::Vector<amrex::DistributionMapping> dm(nlev);
    amrex::Vector<amrex::Geometry> geom(nlev);

    for (int lev = 0; lev < nlev; ++lev) {
        grids[lev] = boxArray(lev);
        dm[lev]    = DistributionMap(lev);
        geom[lev]  = Geom(lev);
    }

    // Rebuild the two sets of masks
    masks = FieldSolver::getLevelMasks(grids, dm, geom);

    gather_masks = FieldSolver::getLevelMasks(
        grids, dm, geom, n_buffer + 1
    );
}

void
GravityPICAmr::Evolve()
{
    //amrex::Real time = 0.0;
    amrex::Real a = a_ini;
    amrex::Print() << "a_ini = " << a_ini;
    amrex::Real kick_factor;
    amrex::Real drift_factor;

    Cosmology cosmo;

    // create file for background data 
    std::ofstream a_log;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        a_log.open("scale_factor.dat", std::ios::out);
        a_log << "# step a\n";
        a_log << std::setprecision(15);
    }
    
    //initialize time at a = 1
    amrex::Real time = 2.0 / (3.0*PhysConst::H_0);
    
    // main PIC loop
    for (int step = 0; step <= max_step; ++step) {
        
        // prepare rhs for poisson solve.
        myPC->DepositCharge(GetVecOfPtrs(rhs));

        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::Real max_rho = rhs[lev]->norm0(0);
            amrex::Print() << "  [step " << step << "] max |rho| at level "
                           << lev << " = " << max_rho << "\n";
        }

        FieldSolver::computePhi(GetVecOfConstPtrs(rhs), GetVecOfPtrs(phi),
                                boxArray(), DistributionMap(), Geom(),
                                GetVecOfConstPtrs(masks), bc);

        FieldSolver::computeE(GetVecOfArrOfPtrs(eField), GetVecOfConstPtrs(phi), Geom());

        myPC->FieldGather(GetVecOfArrOfConstPtrs(eField),
                         GetVecOfConstPtrs(gather_masks));

        if ((particle_output_int > 0) && (step % particle_output_int == 0)) {
            myPC->writeParticles(step);
        }

        if ((plot_int > 0) && (step % plot_int == 0)) {
            WritePlotFile(GetVecOfConstPtrs(rhs), GetVecOfConstPtrs(phi), GetVecOfArrOfConstPtrs(eField), *myPC, geom, step);
            myPC->WriteAsciiFile("particles_ascii_" + std::to_string(step));
        }

        amrex::Real da;
        da = a*std::exp(d_ln_a) - a;
        amrex::Print() << " da is " << da << std::endl;
        
        kick_factor = compute_kick_factor(a, da, cosmo) ;
        drift_factor = compute_drift_factor(a, da, cosmo) ;
     
        // move particles
        myPC->Evolve(GetVecOfArrOfConstPtrs(eField), GetVecOfConstPtrs(rhs), kick_factor, drift_factor);

        // output background data 
        if (amrex::ParallelDescriptor::IOProcessor()) {
            a_log << step << "  " << a << "\n";
            a_log.flush(); 
        }

        amrex::Print() << " max_level " << max_level << std::endl;
        // regrid after evolving and before calling redistribute
        if (max_level > 0 && regrid_int > 0 && step > 0 && step % regrid_int == 0)
        {
            regrid(0, time);
            RebuildMasks();
        }

        if (max_level > 0 && regrid_int > 0 && step > 0 && step % regrid_int == 0)
        {
            regrid(0, time);
            rhs.resize(finest_level + 1);
            phi.resize(finest_level + 1);
            eField.resize(finest_level + 1);
            RebuildMasks();
        }
        
        // Redistribute correctly places particles on new grid
        myPC->Redistribute();

        //try logarithmic timestepping
        a = a * std::exp(d_ln_a);
        //time += dt;
    }
    
    if (amrex::ParallelDescriptor::IOProcessor()) {
        a_log.close();
    }

}

void
GravityPICAmr::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{

    amrex::Print() << "Error EST called !!! " << std::endl;
    const amrex::MultiFab& rho = *rhs[lev];

    for (amrex::MFIter mfi(rho, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox();
        bx.enclosedCells();
        
        const auto& rho_arr = rho.const_array(mfi);
        auto tag_arr = tags.array(mfi);
        const amrex::Real threshold = rho_err;
 
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {   
            const amrex::Real rho_cc =
                    0.125 * (
                        rho_arr(i  , j  , k  ) +
                        rho_arr(i+1, j  , k  ) +
                        rho_arr(i  , j+1, k  ) +
                        rho_arr(i+1, j+1, k  ) +
                        rho_arr(i  , j  , k+1) +
                        rho_arr(i+1, j  , k+1) +
                        rho_arr(i  , j+1, k+1) +
                        rho_arr(i+1, j+1, k+1)
                    );

            //amrex::Print() << " max rho_cc:" << rho_cc << std::endl;

            // using the node at cell's low corner as a stand-in for the whole cell's density
            // averaging the surrounding nodes would be more accurate
            if (std::abs(rho_cc) > threshold) {
                tag_arr(i,j,k) = amrex::TagBox::SET;
            }
        });
    }
}

void
GravityPICAmr::AllocLevelData (int lev, const amrex::BoxArray& ba,
                                const amrex::DistributionMapping& dm)
{
    // Function to allocate a box array of correct regridded size/dimension
    amrex::BoxArray nba = ba;
    nba.surroundingNodes();

    rhs[lev] = std::make_unique<amrex::MultiFab>(nba, dm, Ncomp, 1);
    phi[lev] = std::make_unique<amrex::MultiFab>(nba, dm, Ncomp, 2);
    rhs[lev]->setVal(0.0);
    phi[lev]->setVal(0.0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        eField[lev][idim] = std::make_unique<amrex::MultiFab>(nba, dm, Ncomp, 1);
        eField[lev][idim]->setVal(0.0);
    }
}

// make a new level where there was none
void
GravityPICAmr::MakeNewLevelFromScratch (int lev, amrex::Real /*time*/,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm)
{
    if (static_cast<int>(rhs.size()) <= lev) {
        rhs.resize(lev+1);
        phi.resize(lev+1);
        eField.resize(lev+1);
    }
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    AllocLevelData(lev, ba, dm);
    
}

// make a new level where there was none -- this is the same as from scratch because 
// rhs/phi/eField are fully recomputed from deposited charge every step
void
GravityPICAmr::MakeNewLevelFromCoarse (int lev, amrex::Real /*time*/,
                                        const amrex::BoxArray& ba,
                                        const amrex::DistributionMapping& dm)
{

    if (static_cast<int>(rhs.size()) <= lev) {
        rhs.resize(lev+1);
        phi.resize(lev+1);
        eField.resize(lev+1);
    }
    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    AllocLevelData(lev, ba, dm);

}

// change the shape of a n existing level's BoxArray/DistributionMapping 
void
GravityPICAmr::RemakeLevel (int lev, amrex::Real /*time*/,
                             const amrex::BoxArray& ba,
                             const amrex::DistributionMapping& dm)
{

    SetBoxArray(lev, ba);
    SetDistributionMap(lev, dm);
    AllocLevelData(lev, ba, dm);

}

void
GravityPICAmr::ClearLevel (int lev)
{
    rhs[lev].reset();
    phi[lev].reset();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        eField[lev][idim].reset();
    }
}

