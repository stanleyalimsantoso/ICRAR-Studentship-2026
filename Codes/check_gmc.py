from gmc_ns_plot import plot_gmc_xy
import numpy as np
import matplotlib.pyplot as plt
import os

import os

current_snap = 0

output_dir = r"D:\ICRAR Studentship Data\Run 82\preprocessed"
os.makedirs(output_dir, exist_ok=True)

file_path = os.path.join(output_dir, f"snap_{current_snap:04d}_gas.npz")
out_png  = r"D:\ICRAR Studentship Data\gmc_pic.png"


current_snap = 0
input_path = r"D:\ICRAR Studentship Data\Run 82\raw\tout.dat"
with open(input_path, "r") as file:
    header = file.readline().split()

    num_bodies = int(header[0])
    current_time=float(header[1])

    print(num_bodies)

    max_snap = int(header[2])

    file.readline()
    file.readline()
    file.readline()

    x, y, z, vx, vy, vz, mj, frah, t = [], [], [], [], [], [], [], [], []
    ns_x, ns_y, ns_z, ns_t, ns_m, ns_vx, ns_vy, ns_vz = [], [], [], [], [], [], [], []
    os_m, dm_m = [], []

    for i in range(num_bodies):
        line = file.readline().split()
        if len(line) < 14:
            continue

        iwas = int(float(line[6]))
        if iwas == 2:
            mj.append(float(line[10]))
        elif iwas == 3:
            ns_m.append(float(line[10]))
        elif iwas == 1:
            os_m.append(float(line[10]))
        elif iwas ==0:
            dm_m.append(float(line[10]))

    mj = np.array(mj)
    ns_m = np.array(ns_m)
    os_m = np.array(os_m)
    dm_m = np.array(dm_m)


    np.savez_compressed(
        file_path,
        mj=mj,
        ns_m=ns_m,
        os_m=os_m,
        dm_m=dm_m
    )


    print("Saved snapshot", current_snap, "with", (np.size(mj) + np.size(ns_m) + np.size(os_m) + np.size(dm_m)), "gas particles")




data = np.load(file_path)

m = data["mj"]
ns_m = data["ns_m"]
os_m = data["os_m"]
dm_m = data["dm_m"]

# ns_m = 5.4e-6
# m = 3e-7

m_res = m[0]
os_m_res = os_m[0]
dm_m_res = dm_m[0]

print(np.size(m))
tot_m = np.sum(m)
tot_os_m = np.sum(os_m)
tot_dm_m = np.sum(dm_m)

print(f"stellar mass is: {tot_os_m*6*10**10:e}")
print(f"dark matter mass is {tot_dm_m*6*10**10:e}")
print(f"gas mass is: {tot_m*6*10**10:e}")

print(f"gas mass resolution is: {m_res*6*10**10:e}")
print(f"old star mass resolution is: {os_m_res*6*10**10:e}")
print(f"dark matter star mass resolution is: {dm_m_res*6*10**10:e}")