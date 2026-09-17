import numpy as np 
import scipy
import h5py

rng = np.random.default_rng()
FloatType = np.float32  # double precision: np.float64, for single use np.float32
IntType = np.int32

out_file = "output_two/two_part.txt"

N=2

a_init = 1
Omega_m = 1.0

#### NEW  --- initialize comoving momenta such that rdot = 0 -- true turnaround condition!
#define unit system for momenta
unit_length_in_cm = 3.086e+24
unit_vel_in_cm_s = 1e5
unit_time_in_s = unit_length_in_cm / unit_vel_in_cm_s
       
#convert H_0 to inverse code time by comparing time to reference unit time in inverse hubble units

Mpc_in_cm = 3.086e+24
km_s_in_cm_s = 1e5

ref_unit_time = Mpc_in_cm / km_s_in_cm_s # [Mpc/h/km/s]
H_0 = 100 * ref_unit_time / unit_time_in_s # now H_0 is in 1/[code time]

def Hubble(a):
    # assume pure MD
    return H_0 * np.sqrt(Omega_m/(a*a*a))

pos1 = np.array([[0.5, 1, 1]]) / a_init
pos2 = np.array([1.5, 1, 1]) / a_init

# initialize xdot a^2
vel1 = -a_init**2 * Hubble(a_init)*(pos1 - np.array([1,1,1]))
vel2 = -a_init**2 * Hubble(a_init)*(pos2 - np.array([1,1,1]))

#with respect to COM ...

final_pos = np.vstack([pos1, pos2])
final_vel = np.vstack([vel1, vel2])

weights = np.array([[1.0], [1.0]])

data = np.hstack([final_pos, final_vel, weights])

with open(out_file, "w") as f:
    f.write(f"{N}\n")
    for row in data:
        f.write(" ".join(f"{val:.17g}" for val in row) + "\n")

print(f"Wrote {N} particles to {out_file}")