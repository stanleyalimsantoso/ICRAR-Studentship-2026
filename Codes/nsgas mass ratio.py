import os
import numpy as np
import matplotlib.pyplot as plt

run = 5
M_unit = 6e10
algo = "cdktree"
max_snap = 5
total_gas = False

if total_gas == True:
    tag = "totalgas"
    new_tag = "Total Gas"
else:
    tag = "h2gas"
    new_tag = "H2 Gas"

ns_n = []
gmc_m_wns = []

base_dir = rf"D:\ICRAR Studentship Data\Run {run}"

for desired_snap in range(1,max_snap+1):
    data_dir = os.path.join(base_dir, "postprocessed")
    data_dir = os.path.join(data_dir, f"{algo}_snap_{desired_snap:04d}_{tag}.npz")

    data = np.load(data_dir)
    ns_labels = data["ns_labels"]
    labels = data["labels"]
    ns_x = data["ns_x"]
    gmc_m = data["gmc_m"]
    ns_m = data["ns_m"]

    valid_ns = (ns_labels >= 0)
    if valid_ns.any():
        ns_counts = np.bincount(ns_labels[valid_ns], minlength=gmc_m.size)
        gmcs_with_ns = np.where(ns_counts > 0)[0]
    else:
        gmcs_with_ns = np.array([], dtype=int)

    for gmc_id in gmcs_with_ns:
        mask_ns  = (ns_labels == gmc_id)

        gmc_mass = float(gmc_m[gmc_id])*M_unit

        gmc_m_wns.append(gmc_mass)
        ns_n_num = ns_x[mask_ns]
        ns_m_masked = ns_m[mask_ns]
        ns_m_sum = ns_m_masked.sum()

        ns_n.append(ns_m_sum*M_unit/gmc_mass)

    if gmcs_with_ns.size == 0:
        print(f"Snapshot {desired_snap}: no GMCs with new stars, skipping.")
        continue

gmc_m_wns = np.array(gmc_m_wns)

slope, intercept = np.polyfit(gmc_m_wns, ns_n, 1)
line_y = slope * gmc_m_wns + intercept


plt.plot(gmc_m_wns, line_y, color='red', label=f'Fitted Line (y = {slope:.2g}x + {intercept:.2g})')
plt.scatter(gmc_m_wns, ns_n, s=7, label="data points")
plt.title(f"Ratio of New Star Mass/{new_tag} Mass")
plt.xlabel(f"{new_tag} mass in Solar Mass")
plt.ylabel(f"Ratio of New Star mass and {new_tag} mass")
plt.legend()

output_path = os.path.join(base_dir, "plots", f"Ratio of New star over {new_tag} mass plot.png")
plt.savefig(output_path, dpi=300, bbox_inches="tight")

plt.show()
