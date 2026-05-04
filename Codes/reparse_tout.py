import os
import numpy as np
import argparse

def reparse_tout (input_path, output_path, run):

    r_thr = 0.0114 # 0.0114 for 200pc
    total_gas = False #placeholder only

    x_com = 0
    y_com = 0
    z_com = 0

    base_dir = f"/mnt/d/ICRAR Studentship Data/Run {run}"
    target_path = os.path.join(base_dir, "output", input_path)

    get_data_from_filename = input_path.split("_")
    total_gas_check = get_data_from_filename[0]
    timestep = get_data_from_filename[2].lstrip("0")
    timestep = int(timestep)
    gmc_id = get_data_from_filename[3]

    if total_gas_check == "h2gas":
        total_gas = False
        tag = "h2gas"
    else:
        total_gas = True
        tag = "totalgas"

        

    with open (target_path, "r") as file:
        x_list = []
        y_list = []
        z_list = []
        m_list = []

        header = file.readline().split()
        num_bodies = int(float(header[0]))

        for i in range(num_bodies):
            data = file.readline().split()
            x_list.append(float(data[0]))
            y_list.append(float(data[1]))
            z_list.append(float(data[2]))
            m_list.append(float(data[7]))

        x_list = np.array(x_list)
        y_list = np.array(y_list)
        z_list = np.array(z_list)
        m_list = np.array(m_list)

        x_com = np.sum(x_list * m_list)/sum(m_list)
        y_com = np.sum(y_list * m_list)/sum(m_list)
        z_com = np.sum(z_list * m_list)/sum(m_list)

    npz_path = os.path.join(base_dir, "preprocessed",
                            f"snap_{timestep:04d}_gas.npz")

    data = np.load(npz_path)

    x, y, z, vx, vy, vz, mj, frah, t = [], [], [], [], [], [], [], [], []
    ns_x, ns_y, ns_z, ns_vx, ns_vy, ns_vz, ns_m, ns_t = [], [], [], [], [], [], [], []

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
    ns_m = data["ns_m"]
    ns_t = data["ns_t"]

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    vx = np.array(vx)
    vy = np.array(vy)
    vz = np.array(vz)
    mj = np.array(mj)
    frah = np.array(frah)
    t = np.array(t)

    if total_gas == False:
        mj = mj*frah

    ns_x = np.array(ns_x)
    ns_y = np.array(ns_y)
    ns_z = np.array(ns_z)
    ns_vx = np.array(ns_vx)
    ns_vy = np.array(ns_vy)
    ns_vz = np.array(ns_vz) 
    ns_t = np.array(ns_t)
    ns_m = np.array(ns_m)

    r_squared = (x-x_com)**2 + (y-y_com)**2 + (z-z_com)**2
    r_ns_squared = (ns_x-x_com)**2 + (ns_y-y_com)**2 + (ns_z-z_com)**2

    gasmask = r_squared <= r_thr**2
    starmask = r_ns_squared <= r_thr**2

    gmc_x = np.concatenate((x[gasmask], ns_x[starmask]))
    gmc_y = np.concatenate((y[gasmask], ns_y[starmask]))
    gmc_z = np.concatenate((z[gasmask], ns_z[starmask]))
    gmc_vx = np.concatenate((vx[gasmask], ns_vx[starmask]))
    gmc_vy = np.concatenate((vy[gasmask], ns_vy[starmask]))
    gmc_vz = np.concatenate((vz[gasmask], ns_vz[starmask]))
    gmc_m = np.concatenate((mj[gasmask], ns_m[starmask]))
    gmc_t = np.concatenate((t[gasmask], ns_t[starmask]))

    n_gas  = np.count_nonzero(gasmask)
    n_star = np.count_nonzero(starmask)

    gmc_iwas = np.concatenate((
        2 * np.ones(n_gas, dtype=int),
        3 * np.ones(n_star, dtype=int)
    ))


    outfile = os.path.join(
        output_path,
        f"{tag}_gmc_{timestep:04d}_{gmc_id}"
    )

    with open(outfile, "w") as f:
        f.write(f"{len(gmc_x):d} {float(timestep):.6e}\n")
        for i in range(len(gmc_x)):
            f.write(
                f"{gmc_x[i]:.6e} {gmc_y[i]:.6e} {gmc_z[i]:.6e} "
                f"{gmc_vx[i]:.6e} {gmc_vy[i]:.6e} {gmc_vz[i]:.6e} "
                f"{gmc_iwas[i]:d} {gmc_m[i]:.6e} {gmc_t[i]:.6e}\n"
            )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--run")

    args = p.parse_args()
    reparse_tout(args.input, args.output, args.run)

if __name__ == "__main__":
    main()


        