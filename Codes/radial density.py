import os
import numpy as np
import matplotlib.pyplot as plt
import math

run = 1
desired_snap = 2
algo = "cdktree"

base_dir = f"D:/ICRAR Studentship Data/Run {run}"

gmc_dir = os.path.join(base_dir, "gmc")
gmc_path = os.path.join(gmc_dir, f"{algo}_snap_{desired_snap:04d}_gmcs.npz")

data = np.load(gmc_path)

x = data["x"]
y = data["y"]
z = data["z"]
vx = data["vx"]
vy = data["vy"]
vz = data["vz"]
labels = data["labels"]
gmc_x = data["gmc_x"]
gmc_y = data["gmc_y"]
gmc_z = data["gmc_z"]

gas_dir = os.path.join(base_dir, "preprocessed")
gas_path = os.path.join(gas_dir, f"snap_{desired_snap:04d}_gas.npz")
gas_data = np.load(gas_path)
mj = gas_data["mj"]
frah = gas_data["frah"]
m_all = mj * frah

plots_dir = os.path.join(base_dir, "plots", algo)
os.makedirs(plots_dir, exist_ok=True)

unique_gmcs = np.unique(labels)

for gmc_id in unique_gmcs:
    if gmc_id < 0:
        continue  # skip noise / unassigned if you use -1

    mask = labels == gmc_id

    x_gmc = x[mask]
    y_gmc = y[mask]
    z_gmc = z[mask]
    vx_gmc = vx[mask]
    vy_gmc = vy[mask]
    vz_gmc = vz[mask]
    m_gmc = m_all[mask]

    if m_gmc.size == 0:
        continue

    com_x = gmc_x[gmc_id]
    com_y = gmc_y[gmc_id]
    com_z = gmc_z[gmc_id]

    pos = np.vstack((x_gmc, y_gmc, z_gmc)).T

    dx = pos[:, 0] - com_x
    dy = pos[:, 1] - com_y
    dz = pos[:, 2] - com_z
    r = np.sqrt(dx*dx + dy*dy + dz*dz)

    rad = r.max()
    n_bins = 30
    edges = np.linspace(0.0, rad, n_bins + 1)

    shell_mass, _ = np.histogram(r, bins=edges, weights=m_gmc)

    r_in = edges[:-1]
    r_out = edges[1:]
    shell_vol = (4.0 * math.pi / 3.0) * (r_out**3 - r_in**3)
    shell_mid = 0.5 * (r_in + r_out)

    density = np.where(shell_vol > 0, shell_mass / shell_vol, 0.0)

    plt.figure()
    plt.plot(shell_mid, density, marker="o")
    plt.xlabel("Radius (code units)")
    plt.ylabel("Radial density (mass / volume)")
    plt.title(f"Run {run}, snap {desired_snap}, GMC {gmc_id} ({algo})")
    plt.grid(True)

    out_png = os.path.join(
        plots_dir,
        f"{algo}_snap_{desired_snap:04d}_gmc{int(gmc_id)}.png"
    )

    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    break
