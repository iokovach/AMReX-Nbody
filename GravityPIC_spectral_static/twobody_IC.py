import numpy as np 
import scipy
import h5py

rng = np.random.default_rng()
FloatType = np.float32  # double precision: np.float64, for single use np.float32
IntType = np.int32

out_file = "output/two_part.txt"

N=2

pos1 = np.array([[0.5, 1, 1]])
pos2 = np.array([1.5, 1, 1])

vel1 = np.array([[0,0,0]])
vel2 = np.array([0,0,0])

final_pos = np.vstack([pos1, pos2])
final_vel = np.vstack([vel1, vel2])

weights = np.array([[1.0], [1.0]])

data = np.hstack([final_pos, final_vel, weights])

with open(out_file, "w") as f:
    f.write(f"{N}\n")
    for row in data:
        f.write(" ".join(f"{val:.17g}" for val in row) + "\n")

print(f"Wrote {N} particles to {out_file}")