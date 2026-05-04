import os
import numpy as np
import matplotlib.pyplot as plt

def plot_gmc_xy(file_path, out_png, title=None):
    with open(file_path, "r") as f:
        header = f.readline().strip().split()
        num = int(header[0])
        tnow = float(header[1])

    data = np.loadtxt(file_path, skiprows=1)
    x, y, z, vx, vy, vz, iwas, mass, tt = data.T

    gas_mask = (iwas == 2)
    ns_mask  = (iwas == 3)

    gmc_name = os.path.splitext(os.path.basename(file_path))[0]
    if title is None:
        title = f"{gmc_name} at t = {tnow:.3f}"

    plt.figure(figsize=(6, 4))
    plt.scatter(x[gas_mask], y[gas_mask], s=4, alpha=0.5, label="Gas (iwas=2)")
    plt.scatter(x[ns_mask],  y[ns_mask],  s=30, marker="*", label="New stars (iwas=3)")

    plt.xlabel("x (simulation units)")
    plt.ylabel("y (simulation units)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
