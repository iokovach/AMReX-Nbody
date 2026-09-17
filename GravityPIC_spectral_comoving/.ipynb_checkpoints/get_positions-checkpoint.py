import yt
from yt.frontends import boxlib
from yt.frontends.boxlib.api import AMReXDataset
yt.set_log_level(0)
import numpy as np
import h5py

x_lev0= []
y_lev0 = []
z_lev0 = []

start = 0
stop = 100
step_size = 1

data_dir='/projects/illinois/eng/physics/sheltonj/amrex/GravityPIC/GravityPIC_spectral/output_lna/'

max_x = []
max_y = []
max_z = []
steps_used = []

for step in range(start, stop, step_size):
    print(step)
    ds = AMReXDataset(data_dir + f"plt{step:05d}")

    x = ds.r[("particle0", "particle_position_x")][:]
    y = ds.r[("particle0", "particle_position_y")][:]
    z = ds.r[("particle0", "particle_position_z")][:]  # fixed the y/z bug

    max_x.append(x.max())
    max_y.append(y.max())
    max_z.append(z.max())
    steps_used.append(step)

# stack into one array: shape (n_steps, 4) -> step, max_x, max_y, max_z
out = np.column_stack([steps_used, max_x, max_y, max_z])
np.savetxt("output/max_positions.dat", out, header="step max_x max_y max_z")