import os
import numpy as np
from clustering_algorithms import naive_search, cdktree

r_thr = 0.0029   # 1 length unit is 17.5kpc | 0.0029 for 50pc 0.006 for 100pc
m_thr_upper = 0.00005  # 1 mass unit is 6*10^10 | 0.00017 for 10^7, 0.00005 for 3*10^6
m_thr_lower = 0.000017 # 0.000017 for 10^6

def classify_gmc(gas_dir, gmc_dir, max_timestamp, total_gas, run, algo = cdktree):
    os.makedirs(gmc_dir, exist_ok=True)
    for desired_snap in range(max_timestamp):
        os.makedirs(gmc_dir, exist_ok=True)
        gas_path = os.path.join(gas_dir, f"snap_{desired_snap:04d}_gas.npz")

        if not os.path.exists(gas_path):
            print(f"Run {run}, snapshot {desired_snap}: {gas_path} not found, skipping")
            continue
        data = np.load(gas_path)

        x = data["x"]
        y = data["y"]
        z = data["z"]
        vx = data["vx"]
        vy = data["vy"]
        vz = data["vz"]
        mj = data["mj"]
        frah = data["frah"]
        t = data["t"]
        ns_x = data["ns_x"]
        ns_y = data["ns_y"]
        ns_z = data["ns_z"]
        ns_vx = data["ns_vx"]
        ns_vy = data["ns_vy"]
        ns_vz = data["ns_vz"]
        ns_t = data["ns_t"]
        ns_m = data["ns_m"]

        gmc_x, gmc_y, gmc_z, gmc_vx, gmc_vy, gmc_vz, gmc_m, gmc_n, labels = algo(
            x, y, z, vx, vy, vz, mj, frah, r_thr, m_thr_upper, m_thr_lower, total_gas
        )

        name_algo = algo.__name__

        if total_gas == True:
            out_path = os.path.join(gmc_dir, f"{name_algo}_snap_{desired_snap:04d}_gmcs_totalgas.npz")
        else:
            out_path = os.path.join(gmc_dir, f"{name_algo}_snap_{desired_snap:04d}_gmcs_h2gas.npz")   
        np.savez_compressed(
            out_path,
            x=x,
            y=y,
            z=z,
            vx=vx,
            vy=vy,
            vz=vz,
            mj=mj,
            t=t,
            gmc_x=gmc_x,
            gmc_y=gmc_y,
            gmc_z=gmc_z,
            gmc_vx=gmc_vx,
            gmc_vy=gmc_vy,
            gmc_vz=gmc_vz,
            gmc_m=gmc_m,
            gmc_n=gmc_n,
            labels=labels,
            ns_x=ns_x,
            ns_y=ns_y, 
            ns_z=ns_z,
            ns_vx=ns_vx,
            ns_vy=ns_vy,
            ns_vz=ns_vz,
            ns_t = ns_t,
            ns_m = ns_m
        )

        print(f"Run {run}, snapshot {desired_snap}: {len(gmc_x)} GMCs above mass cut ({name_algo})")
