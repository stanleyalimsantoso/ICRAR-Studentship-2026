import os
import numpy as np
import matplotlib.pyplot as plt
import math
from gmc_ns_plot import plot_gmc_xy

#specific angular momentum -> 2.14 * 10^3 kpc km/s

def parse_name(file_path):
    # expects: tag_gmc_0002_0.dat
    filename = os.path.basename(file_path)
    stem = os.path.splitext(filename)[0]
    parts = stem.split("_")
    tag = parts[0]
    desired_snap = int(parts[2])
    gmc_id = int(parts[3])
    return stem, tag, desired_snap, gmc_id

def read_gmc_dat(file_path):
    with open(file_path, "r") as f:
        header = f.readline().strip().split()
        num = int(header[0])
        tnow = float(header[1])

    data = np.loadtxt(file_path, skiprows=1)
    x, y, z, vx, vy, vz, iwas, mass, tt = data.T
    return num, tnow, x, y, z, vx, vy, vz, iwas, mass, tt

def calculate_star_masses(desired_snap, iwas, mass, tt):

    tnow = desired_snap * 0.2 #very important multiply this by the timestep

    ns_mask = (iwas == 3)
    if np.sum(ns_mask) == 0:
        return 0.0, 0.0

    age = tnow - tt[ns_mask]
    m = mass[ns_mask]

    new_mask = (age < 0.71)                 
    agb_mask = (0.14 < age) & (age < 0.71)  

    new_star_mass = float(m[new_mask].sum())
    agb_star_mass = float(m[agb_mask].sum())
    return new_star_mass, agb_star_mass

def plot_radial_density_gas(x, y, z, iwas, mass, out_png, title):
    gas_mask = (iwas == 2)
    x_g = x[gas_mask]; y_g = y[gas_mask]; z_g = z[gas_mask]; m_g = mass[gas_mask]

    if m_g.size == 0:
        return None, None

    m_tot = m_g.sum()
    com_x = (x_g * m_g).sum() / m_tot
    com_y = (y_g * m_g).sum() / m_tot
    com_z = (z_g * m_g).sum() / m_tot

    r = np.sqrt((x_g - com_x)**2 + (y_g - com_y)**2 + (z_g - com_z)**2)
    rad = float(r.max())
    if rad <= 0.0:
        rad = 1.0e-6

    n_bins = 30
    edges = np.linspace(0.0, rad, n_bins + 1)
    shell_mass, _ = np.histogram(r, bins=edges, weights=m_g)

    r_in = edges[:-1]
    r_out = edges[1:]
    shell_vol = (4.0 * math.pi / 3.0) * (r_out**3 - r_in**3)
    shell_mid = 0.5 * (r_in + r_out)

    density = np.where(shell_vol > 0, shell_mass / shell_vol, 0.0)

    mask = (shell_mid > 0) & (density > 0)

    plt.figure()
    plt.plot(shell_mid[mask], density[mask], marker="o")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radius (code units)")
    plt.ylabel("Radial density (mass / volume)")
    plt.title(title)
    plt.grid(True, which="both")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    return shell_mid, density

def specific_angular_momentum_gas(x, y, z, vx, vy, vz, iwas, mass):
    gas_mask = (iwas == 2)
    x_g = x[gas_mask]; y_g = y[gas_mask]; z_g = z[gas_mask]
    vx_g = vx[gas_mask]; vy_g = vy[gas_mask]; vz_g = vz[gas_mask]
    m_g = mass[gas_mask]

    if m_g.size == 0:
        return np.array([0.0, 0.0, 0.0]), 0.0

    m_tot = m_g.sum()

    com_x = (x_g * m_g).sum() / m_tot
    com_y = (y_g * m_g).sum() / m_tot
    com_z = (z_g * m_g).sum() / m_tot

    com_vx = (vx_g * m_g).sum() / m_tot
    com_vy = (vy_g * m_g).sum() / m_tot
    com_vz = (vz_g * m_g).sum() / m_tot

    rx = x_g - com_x
    ry = y_g - com_y
    rz = z_g - com_z

    vx_rel = vx_g - com_vx
    vy_rel = vy_g - com_vy
    vz_rel = vz_g - com_vz

    Lx = (m_g * (ry * vz_rel - rz * vy_rel)).sum()
    Ly = (m_g * (rz * vx_rel - rx * vz_rel)).sum()
    Lz = (m_g * (rx * vy_rel - ry * vx_rel)).sum()

    j_vec = np.array([Lx, Ly, Lz]) / m_tot
    j_mag = float(np.sqrt((j_vec * j_vec).sum()))
    return j_vec, j_mag

def build_profile(file_path, profiles_root):
    gmc_name, tag, desired_snap, gmc_id = parse_name(file_path)

    num, tnow, x, y, z, vx, vy, vz, iwas, mass, tt = read_gmc_dat(file_path)
    gas_mask = (iwas == 2)

    gmc_dir = os.path.join(profiles_root, gmc_name)
    os.makedirs(gmc_dir, exist_ok=True)

    radial_png = os.path.join(gmc_dir, "radial_density.png")
    xy_png = os.path.join(gmc_dir, "gmc_xy.png")
    summary_path = os.path.join(gmc_dir, "profile.txt")

    new_star_mass, agb_star_mass = calculate_star_masses(desired_snap, iwas, mass, tt)
    total_gas_mass = float(mass[gas_mask].sum())
    j_vec, j_mag = specific_angular_momentum_gas(x, y, z, vx, vy, vz, iwas, mass)

    if total_gas_mass > 0:
        agb_frac = agb_star_mass/total_gas_mass
    else:
        agb_frac = 0

    plot_radial_density_gas(
        x, y, z, iwas, mass, radial_png,
        title=f"{gmc_name} radial density"
    )
    plot_gmc_xy(file_path, xy_png)

    with open(summary_path, "w") as f:
        f.write(f"gmc_name {gmc_name}\n")
        f.write(f"tnow {tnow:.6e}\n")
        f.write(f"total_gas_mass {total_gas_mass:.6e}\n")
        f.write(f"new_star_mass(<0.1Gyr) {new_star_mass:.6e}\n")
        f.write(f"agb_star_mass(0.02-0.1Gyr) {agb_star_mass:.6e}\n")
        f.write(f"fraction of agb mass to gas mass {agb_frac:.6e}\n")
        f.write(f"specific_angular_momentum_jx {j_vec[0]:.6e}\n")
        f.write(f"specific_angular_momentum_jy {j_vec[1]:.6e}\n")
        f.write(f"specific_angular_momentum_jz {j_vec[2]:.6e}\n")
        f.write(f"specific_angular_momentum_mag {j_mag:.6e}\n")

    return agb_frac