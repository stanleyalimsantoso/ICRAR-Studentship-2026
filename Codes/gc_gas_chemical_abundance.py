import os
import numpy as np
import matplotlib.pyplot as plt

max_timestamp = 10
run = 72

base_path = f"D:\ICRAR Studentship Data\Run {run}\GC Data"
preprocessed_dir = os.path.join(base_path, "preprocessed")
plots_dir = os.path.join(base_path, "plots")

for desired_snap in range(max_timestamp):

    input_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")
    data = np.load(input_path)

    gas_K = data["gas_K"]
    gas_Mg = data["gas_Mg"]
    gas_Fe = data["gas_Fe"]

    Mg_Fe_gas = (np.log10(gas_Mg/gas_Fe) - np.log10(5.15e-4/1.17e-3))
    K_Fe_gas  = (np.log10(gas_K/gas_Fe)  - np.log10(3.649e-06/1.17e-3))

    output_path = os.path.join(plots_dir, f"chemical abundance snap {desired_snap}")

    plt.close("all")
    plt.figure()
    plt.scatter(Mg_Fe_gas, K_Fe_gas, s=2, c="green", label="gas")
    plt.xlabel("[Mg/Fe]")
    plt.ylabel("[K/Fe]")
    plt.ylim(-0.2, 1.0)
    plt.xlim(-0.5, 0.7)
    plt.title(f"Mg/Fe vs K/Fe (snap {desired_snap:04d})")
    plt.tight_layout()
    plt.legend()
    plt.savefig(output_path, dpi=200)
    plt.close()


