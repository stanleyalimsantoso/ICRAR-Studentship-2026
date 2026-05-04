import matplotlib.pyplot as plt
import numpy as np
import os


def heatmap(gmc_dir, plots_dir, run, max_timestamp, total_gas, algo="cdktree"):

    for snap in range(max_timestamp):
        if total_gas:
            gmc_path = os.path.join(gmc_dir, f"{algo}_snap_{snap:04d}_gmcs_totalgas.npz")
            tag = "totalgas"
        else:
            gmc_path = os.path.join(gmc_dir, f"{algo}_snap_{snap:04d}_gmcs_h2gas.npz")
            tag = "h2gas"

        if not os.path.exists(gmc_path):
            print(f"Run {run}, snapshot {snap}: {gmc_path} not found, skipping")
            continue

        data = np.load(gmc_path)

        gmc_x_com = data["gmc_x"]
        gmc_y_com = data["gmc_y"]

        x = data["x"]
        y = data["y"]

        bins = 200
        h = plt.hist2d(x, y, bins, cmap="plasma", density=False)
        plt.scatter(gmc_x_com, gmc_y_com, s=5, c="white")
        plt.colorbar(h[3])

        os.makedirs(plots_dir, exist_ok=True)

        out_png = os.path.join(plots_dir, f"{algo}_snap_{snap:04d}_{tag}_wgmc.png")

        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Created heatmap of timestep: {snap} with detected {len(gmc_x_com)} GMCs")

    for snap in range(max_timestamp+1):

        if total_gas:
            gmc_path = os.path.join(gmc_dir, f"{algo}_snap_{snap:04d}_gmcs_totalgas.npz")
            tag = "totalgas"
        else:
            gmc_path = os.path.join(gmc_dir, f"{algo}_snap_{snap:04d}_gmcs_h2gas.npz")
            tag = "h2gas"

        if not os.path.exists(gmc_path):
            print(f"Run {run}, snapshot {snap}: {gmc_path} not found, skipping")
            continue

        data = np.load(gmc_path)

        gmc_x_com = data["gmc_x"]
        gmc_y_com = data["gmc_y"]

        x = data["x"]
        y = data["y"]

        bins = 200
        h = plt.hist2d(x, y, bins, cmap="plasma", density=False)
        plt.colorbar(h[3])

        os.makedirs(plots_dir, exist_ok=True)

        out_png = os.path.join(plots_dir, f"{algo}_snap_{snap:04d}_{tag}.png")

        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Created heatmap of timestep: {snap}")