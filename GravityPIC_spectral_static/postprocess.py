import numpy as np
import h5py
import glob
import sys
from scipy.interpolate import CubicSpline, interp1d
import re
import os

#usage python postprocess.py [path to output dir]
data_path=sys.argv[1]

start = 0
stop = 1000

with open(os.path.join(data_path, "inputs")) as f:
    dt = float(next(line.split("=")[1] for line in f if line.startswith("dt =")))

print("found dt =", dt)

pos=[]
vel=[]
acc=[]
ids=[]
mass=[]

def snap_num(path):
    return int(re.search(r"particles_ascii_(\d+)", path).group(1))

files = sorted(glob.glob(data_path + r"particles_ascii_*"), key=snap_num)

steps = np.array([snap_num(f) for f in files])
times = steps * dt


for f in files:
    data = np.loadtxt(f, skiprows=5)
    ids.append(data[:, 3].astype(np.int64))
    mass = data[:, 5]
    pos.append(data[:, 0:3])   
    vel.append(data[:, 6:9])
    acc.append(data[:, 9:12])


# find which particles exist at every snapshot
common_ids = set(ids[0])
for id_arr in ids[1:]:
    common_ids &= set(id_arr)
common_ids = np.sort(np.array(list(common_ids)))
print(f"{len(common_ids)} / {len(ids[0])} particles survive all {len(files)} snapshots")

pos_aligned, vel_aligned, acc_aligned = [], [], []

for snap in range(len(files)):
    order = np.argsort(ids[snap])
    sel = order[np.searchsorted(ids[snap][order], common_ids)]
    pos_aligned.append(pos[snap][sel])
    vel_aligned.append(vel[snap][sel])
    acc_aligned.append(acc[snap][sel])
    
pos, vel, acc = pos_aligned, vel_aligned, acc_aligned

#track outermost particle position in x/y/z
outfile = h5py.File(data_path+'outermost_part.h5', 'w')
outfile.create_dataset("time", data=times)
outfile.create_dataset("x", data= np.max(np.array(pos)[...,0], axis=0))
outfile.create_dataset("y", data= np.max(np.array(pos)[...,1], axis=0))
outfile.create_dataset("z", data= np.max(np.array(pos)[...,2], axis=0))

#compute GW quantities

Qs=[]
Qddots=[]
Qdddots=[]

print("mass [0] shape", mass[0].shape)
def quadrupole_ddot(y, ydot, yddot, mass):
    outer = (np.einsum('ni,nj->ij', yddot, y)
             + np.einsum('ni,nj->ij', y, yddot)
             + 2 * np.einsum('ni,nj->ij', ydot, ydot))
    
    Qij_ddot = mass[0] * outer

    trace_term = mass[0] * (2/3) * np.sum(ydot**2 + y * yddot)
    Qij_ddot -= trace_term * np.eye(3)
    
    return Qij_ddot

def quadrupole_dddot(y, ydot, yddot, ydddot, mass):
    outer = (np.einsum('ni,nj->ij', y, ydddot)
             + np.einsum('ni,nj->ij', ydddot, y)
             + 3*np.einsum('ni,nj->ij', yddot, ydot)
             + 3*np.einsum('ni,nj->ij', ydot, yddot))
    
    Qij_dddot = mass[0] * outer

    trace_term = mass[0] * (2/3) * (np.sum(y * ydddot + 3*ydot*yddot))
    Qij_dddot -= trace_term * np.eye(3)
    
    return Qij_dddot

#now we need to interpolate for the third derivative:
#xdddot=CubicSpline(times, np.array(acc), axis=0).derivative()
#cubic spline was expensive in memory so let's try gradient 
acc_array = np.stack(acc, axis=0)
xdddot_array = np.gradient(acc_array, times, axis=0, edge_order=2)

for snap in np.arange(0,len(times)):
    y = pos[snap] - np.average(pos[snap], axis=0)

    N = y.shape[0]
    outer_sum = np.einsum('ni,nj->ij', y, y)
    r2_sum = np.sum(y**2)                        
    
    Qij = outer_sum - (1/3) * r2_sum * np.eye(3)

    Qs.append(Qij)
    
    #get r/v/a for pcls rel to COM
    y=pos[snap] - np.average(pos[snap], axis=0)
    ydot=vel[snap] - np.average(vel[snap], axis=0)
    yddot=acc[snap] - np.average(acc[snap], axis=0)
    ydddot=xdddot_array[snap] - np.average(xdddot_array[snap], axis=0)

    Qddots.append(quadrupole_ddot(y, ydot, yddot, mass))
    Qdddots.append(quadrupole_dddot(y, ydot, yddot, ydddot, mass))

Q_final=np.array(Qs)
Qddot_final=np.array(Qddots)
Qdddot_final=np.array(Qdddots)

##FFT
newtimes= np.linspace(times[0],times[-1],times.shape[0])
dt=newtimes[1]-newtimes[0]

#interpolate to a regular time grid
Q_interp=interp1d(times, Q_final, axis=0)

# real FFT since signal is real
Q_f = np.fft.rfft(Q_interp(newtimes), axis=0)
freqs = np.fft.rfftfreq(len(newtimes), d=dt)

#interpolate to a regular time grid
Qddot_interp=interp1d(times, Qddot_final, axis=0)

# real FFT since signal is real
Qddot_f = np.fft.rfft(Qddot_interp(newtimes), axis=0)
freqs = np.fft.rfftfreq(len(newtimes), d=dt)

#interpolate to a regular time grid
Qdddot_interp=interp1d(times, Qdddot_final, axis=0)

# real FFT since signal is real
Qdddot_f = np.fft.rfft(Qdddot_interp(newtimes), axis=0)
freqs = np.fft.rfftfreq(len(newtimes), d=dt)

power = freqs*np.sum(np.abs(Qdddot_f)**2, axis=(1,2))

outfile = h5py.File(data_path+'GW_diagnostics.h5', 'w')
outfile.create_dataset("times", data=times)
outfile.create_dataset("Q", data=Q_final)
outfile.create_dataset("Qddot", data=Qddot_final)
outfile.create_dataset("Qdddot", data=Qdddot_final)

outfile.create_dataset("freq", data=freqs)
outfile.create_dataset("QFFT", data=Q_f)
outfile.create_dataset("QddFFT", data=Qddot_f)
outfile.create_dataset("QdddFFT", data=Qdddot_f)
outfile.create_dataset("power", data=power)

#we can also compute half mass radius

Nsnap = len(pos)
r50 = np.zeros(Nsnap)
#halfmass radius
for snap in range(Nsnap):
    y = pos[snap] - np.average(pos[snap], axis=0)
    # Half-mass radius
    r_i = np.linalg.norm(y, axis=1)
    order = np.argsort(r_i)
    m_sorted = mass[order]
    cum_mass = np.cumsum(m_sorted)
    idx50 = np.searchsorted(cum_mass, 0.5*cum_mass[-1])
    r50[snap] = r_i[order][idx50]

np.savetxt(data_path+'r50.txt', r50)






