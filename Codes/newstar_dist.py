import os
import numpy as np

t_thr_upper = 1.4 #1 unit is 0.14Gyr | 0.7 for 0.01gyr
t_thr_lower = 0.14 # 0.14 for 0.02gyr
r_thr_upper = 0.0114 #1 unit is 17.5 kpc | 0.006 for 100pc 0.0114 for 200pc
r_thr_lower = 0.0029

def newstar_dist(target_dir, output_dir, max_timestamp, run, total_gas, algo="cdktree"):


    max_snap = max_timestamp

    for desired_snap in range(max_snap+1):
        if total_gas:
            specific_target_dir = os.path.join(target_dir, f"{algo}_snap_{desired_snap:04d}_gmcs_totalgas.npz")
        else:
            specific_target_dir = os.path.join(target_dir, f"{algo}_snap_{desired_snap:04d}_gmcs_h2gas.npz")

        gmc_with_agb = 0

        if not os.path.exists(specific_target_dir):
            print(f"Snapshot {desired_snap}: {specific_target_dir} not found, skipping.")
            continue

        data = np.load(specific_target_dir)
        x = data["x"]
        y = data["y"]
        z = data["z"]
        vx = data["vx"]
        vy = data["vy"]
        vz = data["vz"]
        mj = data["mj"]
        t = data["t"]
        ns_x = data["ns_x"]
        ns_y = data["ns_y"]
        ns_z = data["ns_z"]
        ns_t = data["ns_t"]
        ns_m = data["ns_m"]
        ns_vx = data["ns_vx"]
        ns_vy = data["ns_vy"]
        ns_vz = data["ns_vz"]
        gmc_x = data["gmc_x"]
        gmc_y = data["gmc_y"]
        gmc_z = data["gmc_z"]
        gmc_vx = data["gmc_vx"]
        gmc_vy = data["gmc_vy"]
        gmc_vz = data["gmc_vz"]
        gmc_m = data["gmc_m"]
        labels = data["labels"]


        if gmc_x.size == 0:
            print(f"Run {run}, snapshot {desired_snap}: gmc_x is empty (0 GMCs). Skipping.")
            continue

        nstar = len(ns_x)
        ns_age = desired_snap*0.2 - ns_t
        ns_labels = -1 * np.ones(nstar, dtype=int)

        unique_gmcs = np.unique(labels)
        unique_gmcs = unique_gmcs[unique_gmcs != -1]

        unique_gmcs = unique_gmcs[(unique_gmcs >= 0) & (unique_gmcs < gmc_x.size)]

        nstar = len(ns_x)
        visited = np.zeros(nstar, dtype=bool)
        total_m = 0.0

        for gmc_id in unique_gmcs:
            com_x = gmc_x[gmc_id]
            com_y = gmc_y[gmc_id]
            com_z = gmc_z[gmc_id]

            dx = ns_x - com_x
            dy = ns_y - com_y
            dz = ns_z - com_z
            dist = np.sqrt(dx*dx + dy*dy + dz*dz)

            mask_too_young = (dist <= r_thr_lower) & (ns_age < t_thr_lower)
            if np.any(mask_too_young):
                continue

            mask = (~visited) & (dist <= r_thr_upper) & (ns_age <= t_thr_upper) & (ns_age >= t_thr_lower)
            idx = np.where(mask)[0]

            if idx.size == 0:
                continue

            gmc_with_agb += 1

            ns_labels[idx] = gmc_id

            visited[idx] = True
            total_m += ns_m[idx].sum()

        os.makedirs(output_dir, exist_ok=True)

        if total_gas == True:
            output_specific = os.path.join(output_dir, f"{algo}_snap_{desired_snap:04d}_totalgas.npz")
        else:
            output_specific = os.path.join(output_dir, f"{algo}_snap_{desired_snap:04d}_h2gas.npz")

        np.savez_compressed(
            output_specific,
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
            labels=labels,
            ns_x=ns_x,
            ns_y=ns_y, 
            ns_z=ns_z,
            ns_vx=ns_vx,
            ns_vy=ns_vy,
            ns_vz=ns_vz,
            ns_t = ns_t,
            ns_m = ns_m,
            nstar = nstar,
            ns_labels = ns_labels
        )

        n_selected = int(np.sum(visited))
        print(
        f"Run {run}, snapshot {desired_snap}: "
            f"{gmc_with_agb} GMCs with AGB stars; "
            f"{n_selected} AGB stars; "
            f"total AGB star mass = {total_m:.6e}")