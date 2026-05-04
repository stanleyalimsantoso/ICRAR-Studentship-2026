import os
import numpy as np
import matplotlib.pyplot as plt

max_timestamp = 5
run = 18

base_path = f"D:\ICRAR Studentship Data\Run {run}\GC Data"
preprocessed_dir = os.path.join(base_path, "preprocessed")
radial_dir = os.path.join(base_path, "radial")


def specific_angular_momentum(x,y,z,vx,vy,vz,m):
    px, py, pz = m*vx, m*vy, m*vz
    r = np.column_stack((x, y, z))
    p = np.column_stack((px, py, pz))

    L = np.cross(r, p)
    Ltot = np.sum(L, axis=0)

    Lmag = np.sqrt(Ltot[0]**2 + Ltot[1]**2 + Ltot[2]**2)
    mtot = np.sum(m)

    return Lmag / mtot
for desired_snap in range(max_timestamp):
    snap_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")

    data = np.load(snap_path)

    x = data["x"]; y = data["y"]; z = data["z"]; vx = data["vx"]; vy = data["vy"]; vz = data["vz"]; mj = data["mj"]
    ns_x = data["ns_x"]; ns_y = data["ns_y"]; ns_z = data["ns_z"]; ns_vx = data["ns_vx"]; ns_vy = data["ns_vy"]; ns_vz = data["ns_vz"]; ns_m = data["ns_m"]
    os_x = data["os_x"]; os_y = data["os_y"]; os_z = data["os_z"]; os_vx = data["os_vx"]; os_vy = data["os_vy"]; os_vz = data["os_vz"]; os_m = data["os_m"]


    total_m = np.sum(mj) + np.sum(ns_m) + np.sum(os_m)

    com_x = (np.sum(x*mj) + np.sum(ns_m*ns_x) + np.sum(os_m*os_x)) / total_m
    com_y = (np.sum(y*mj) + np.sum(ns_m*ns_y) + np.sum(os_m*os_y)) / total_m
    com_z = (np.sum(z*mj) + np.sum(ns_m*ns_z) + np.sum(os_m*os_z)) / total_m

    x = x - com_x; y = y - com_y; z = z - com_z
    ns_x = ns_x - com_x; ns_y = ns_y - com_y; ns_z = ns_z - com_z
    os_x = os_x - com_x; os_y = os_y - com_y; os_z = os_z - com_z

    com_vx = (np.sum(vx*mj) + np.sum(ns_m*ns_vx) + np.sum(os_m*os_vx)) / total_m
    com_vy = (np.sum(vy*mj) + np.sum(ns_m*ns_vy) + np.sum(os_m*os_vy)) / total_m
    com_vz = (np.sum(vz*mj) + np.sum(ns_m*ns_vz) + np.sum(os_m*os_vz)) / total_m

    vx = vx - com_vx; vy = vy - com_vy; vz = vz - com_vz
    ns_vx = ns_vx - com_vx; ns_vy = ns_vy - com_vy; ns_vz = ns_vz - com_vz
    os_vx = os_vx - com_vx; os_vy = os_vy - com_vy; os_vz = os_vz - com_vz

    gas_sam = specific_angular_momentum(x,y,z,vx,vy,vz,mj)
    ns_sam = specific_angular_momentum(ns_x, ns_y, ns_z, ns_vx, ns_vy, ns_vz, ns_m)
    os_sam = specific_angular_momentum(os_x, os_y, os_z, os_vx, os_vy, os_vz, os_m)

    gas_mean = np.mean(gas_sam) if gas_sam.size > 0 else np.nan
    ns_mean  = np.mean(ns_sam)  if ns_sam.size  > 0 else np.nan
    os_mean  = np.mean(os_sam)  if os_sam.size  > 0 else np.nan

    out_txt = os.path.join(radial_dir, f"mean_specific_angular_momentum_snap_{desired_snap:04d}.txt")

    vals = np.array([[gas_mean, ns_mean, os_mean]], dtype=float)

    np.savetxt(
        out_txt,
        vals,
        fmt="%.8e",
        header="gas_mean_sam  new_stars_mean_sam  old_stars_mean_sam"
    )