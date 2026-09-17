import numpy as np 
from scipy.spatial import cKDTree
import scipy
import h5py

rng = np.random.default_rng()
FloatType = np.float32  # double precision: np.float64, for single use np.float32
IntType = np.int32

R = 1
N = 3000
Mtot = 1
out_file = "output/sphere3000_ics_center1.txt"

#position of center of box
centerx = 1
centery = 1
centerz = 1

center_pos = np.array([centerx,centery,centerz])

Vel = np.zeros((N,3), dtype=FloatType)
ids = np.arange(N)
mass_weight = Mtot / N

space=((R**3)*(4*np.pi/3)/N)**(1/3)

l=[]

while len(l) < N:
    x, y, z = rng.uniform(-R, R, size=3)
    r2 = x*x + y*y + z*z
    if r2 < R**2:
        l.append([x, y, z])
        
final_pos=np.array(l)
final_vel=Vel

weights = np.full((N, 1), mass_weight)

# sphere centered at 0. We want it to be at the center of the box
final_pos += center_pos

data = np.hstack([final_pos, final_vel, weights])

with open(out_file, "w") as f:
    f.write(f"{N}\n")
    for row in data:
        f.write(" ".join(f"{val:.17g}" for val in row) + "\n")

print(f"Wrote {N} particles to {out_file}")