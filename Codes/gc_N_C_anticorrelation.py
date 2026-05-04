import os
import numpy as np
import matplotlib.pyplot as plt

max_timestamp = 10
run = 69

base_path = f"D:\ICRAR Studentship Data\Run {run}\GC Data"
preprocessed_dir = os.path.join(base_path, "preprocessed")
plots_dir = os.path.join(base_path, "plots")

for desired_snap in range(max_timestamp):

    input_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")
    data = np.load(input_path)

    ns_N = data["ns_N"]
    ns_C = data["ns_C"]
    ns_Fe = data["ns_Fe"]

    N_Fe_ns = (np.log10(ns_N/ns_Fe) - np.log10(1.11e-3/1.17e-3))
    C_Fe_ns  = (np.log10(ns_C/ns_Fe)  - np.log10(3.03e-03/1.17e-3))

    output_path = os.path.join(plots_dir, f"chemical abundance snap {desired_snap}")

    plt.close("all")
    plt.figure()
    plt.scatter(N_Fe_ns, C_Fe_ns, s=2, c="red", label="new stars")
    plt.xlabel("[Mg/Fe]")
    plt.ylabel("[K/Fe]")
    plt.ylim(-0.2, 1.0)
    plt.xlim(-0.5, 0.7)
    plt.title(f"Mg/Fe vs K/Fe (snap {desired_snap:04d})")
    plt.tight_layout()
    plt.legend()
    plt.savefig(output_path, dpi=200)
    plt.close()


