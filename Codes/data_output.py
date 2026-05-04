import os
import numpy as np

num_gmcs = 3

def data_output(open_dir, out_dir, max_timestamp, run, total_gas, algo="cdktree"):
    if total_gas:
        tag = "totalgas"
    else:
        tag = "h2gas"

    os.makedirs(out_dir, exist_ok=True)

    for desired_snap in range(max_timestamp):

        specific_open_dir = os.path.join(open_dir,f"{algo}_snap_{desired_snap:04d}_{tag}.npz")

        if not os.path.exists(specific_open_dir):
            print(f"Snapshot {desired_snap}: {specific_open_dir} not found, skipping.")
            continue

        data = np.load(specific_open_dir)

        x_gas_all  = data["x"]
        y_gas_all  = data["y"]
        z_gas_all  = data["z"]
        vx_gas_all = data["vx"]
        vy_gas_all = data["vy"]
        vz_gas_all = data["vz"]
        mj_gas_all = data["mj"]
        t_gas_all  = data["t"]
        labels = data["labels"]

        ns_x_all  = data["ns_x"]
        ns_y_all  = data["ns_y"]
        ns_z_all  = data["ns_z"]
        ns_t_all  = data["ns_t"]
        ns_m_all  = data["ns_m"]
        ns_vx_all = data["ns_vx"]
        ns_vy_all = data["ns_vy"]
        ns_vz_all = data["ns_vz"]
        ns_labels = data["ns_labels"]

        gmc_x  = data["gmc_x"]
        gmc_y  = data["gmc_y"]
        gmc_z  = data["gmc_z"]
        gmc_vx = data["gmc_vx"]
        gmc_vy = data["gmc_vy"]
        gmc_vz = data["gmc_vz"]
        gmc_m  = data["gmc_m"]

        valid_ns = (ns_labels >= 0)
        if valid_ns.any():
            ns_counts = np.bincount(ns_labels[valid_ns], minlength=gmc_m.size)
            gmcs_with_ns = np.where(ns_counts > 0)[0]
        else:
            gmcs_with_ns = np.array([], dtype=int)

        if gmcs_with_ns.size == 0:
            print(f"Snapshot {desired_snap}: no GMCs with new stars, skipping.")
            continue

        gmc_density = np.zeros_like(gmc_m)

        for gmc_id in gmcs_with_ns:
            mask_gas = (labels == gmc_id)
            if not np.any(mask_gas):
                continue

            dx = x_gas_all[mask_gas] - gmc_x[gmc_id]
            dy = y_gas_all[mask_gas] - gmc_y[gmc_id]
            dz = z_gas_all[mask_gas] - gmc_z[gmc_id]

            r = np.sqrt(dx*dx + dy*dy + dz*dz)
            if r.size == 0:
                continue

            mask_core = (r <= 0.00286)
            if not np.any(mask_core):
                continue

            core_mass = mj_gas_all[mask_gas][mask_core].sum()
            core_volume = 9.77e-8 #simulation unit

            gmc_density[gmc_id] = core_mass / core_volume

        sorted_idx = np.argsort(gmc_density[gmcs_with_ns])
        top_gmc_idx = gmcs_with_ns[sorted_idx[-num_gmcs:]]

        for gmc_id in top_gmc_idx:
            mask_gas = (labels == gmc_id)
            mask_ns  = (ns_labels == gmc_id)

            x  = x_gas_all[mask_gas]
            y  = y_gas_all[mask_gas]
            z  = z_gas_all[mask_gas]
            vx = vx_gas_all[mask_gas]
            vy = vy_gas_all[mask_gas]
            vz = vz_gas_all[mask_gas]
            mj = mj_gas_all[mask_gas]
            t  = t_gas_all[mask_gas]

            ns_x  = ns_x_all[mask_ns]
            ns_y  = ns_y_all[mask_ns]
            ns_z  = ns_z_all[mask_ns]
            ns_vx = ns_vx_all[mask_ns]
            ns_vy = ns_vy_all[mask_ns]
            ns_vz = ns_vz_all[mask_ns]
            ns_m  = ns_m_all[mask_ns]
            ns_t  = ns_t_all[mask_ns]

            if ns_x.size == 0:
                continue

            x_all = np.concatenate([x,  ns_x])  - gmc_x[gmc_id]
            y_all = np.concatenate([y,  ns_y])  - gmc_y[gmc_id]
            z_all = np.concatenate([z,  ns_z])  - gmc_z[gmc_id]

            vx_all = np.concatenate([vx, ns_vx]) - gmc_vx[gmc_id]
            vy_all = np.concatenate([vy, ns_vy]) - gmc_vy[gmc_id]
            vz_all = np.concatenate([vz, ns_vz]) - gmc_vz[gmc_id]

            mass_all = np.concatenate([mj, ns_m])
            tt_all   = np.concatenate([t,  ns_t])

            iwas_gas = 2 * np.ones_like(mj,   dtype=int)
            iwas_ns  = 3 * np.ones_like(ns_m, dtype=int)
            iwas_all = np.concatenate([iwas_gas, iwas_ns])

            num  = len(iwas_all)
            tnow = float(desired_snap) 

            outfile = os.path.join(
                out_dir,
                f"{tag}_gmc_{desired_snap:04d}_{gmc_id}.dat"
            )

            print("Writing:", outfile, "num =", num)
            with open(outfile, "w") as f:
                f.write(f"{num:d} {tnow:.6e}\n")
                for i in range(num):
                    f.write(
                        f"{x_all[i]:.6e} {y_all[i]:.6e} {z_all[i]:.6e} "
                        f"{vx_all[i]:.6e} {vy_all[i]:.6e} {vz_all[i]:.6e} "
                        f"{iwas_all[i]:d} {mass_all[i]:.6e} {tt_all[i]:.6e}\n"
                    )
