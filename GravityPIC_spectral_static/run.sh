#!/bin/bash -l
##SBATCH --mail-type=BEGIN,END,FAIL
##SBATCH --mail-user=iok2@illinois.edu
#SBATCH --time=01:00:00
#SBATCH --partition=IllinoisComputes
#SBATCH --account=sheltonj-ic
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --job-name Run_Amrex
#SBATCH --output="logs/run.out" 
#SBATCH --error="logs/run.err" 

ROOT=/projects/illinois/eng/physics/sheltonj/amrex/GravityPIC/GravityPIC_spectral_static/output_1e7Msun_lev1_rho10
cd $ROOT

mkdir logs

/usr/bin/time -v -o logs/benchmark.txt ../main3d.gnu.ex inputs