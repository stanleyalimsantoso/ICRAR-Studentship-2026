import os
import numpy as np
import matplotlib.pyplot as plt

max_timestamps = 10
run=60

for desired_snap in range(max_timestamps):

    base_dir = "D:\ICRAR Studentship Data"
    file_path = os.path.join(base_dir, f"Run {run}", "GC data", "preprocessed", f"snap_{desired_snap:04d}_gas.npz")
    out_png = os.path.join(base_dir, f"Run {run}", "GC data", "plots", f"stars plot {desired_snap:04d}")

    data = np.load(file_path)

    ns_x = data["ns_x"]
    ns_y = data["ns_y"]
    os_x = data["os_x"]
    os_y = data["os_y"]

    ns_mask = ((ns_x > -20) & (ns_x < 20) &
            (ns_y > -20) & (ns_y < 20))

    os_mask = ((os_x > -20) & (os_x < 20) &
            (os_y > -20) & (os_y < 20))


    plt.scatter(ns_x[ns_mask], ns_y[ns_mask], s=4, label="new stars")
    plt.scatter(os_x[os_mask], os_y[os_mask], s=4, label="old stars")

    plt.xlabel("x (simulation units)")
    plt.ylabel("y (simulation units)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()
