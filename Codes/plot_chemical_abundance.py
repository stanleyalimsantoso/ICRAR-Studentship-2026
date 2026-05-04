import matplotlib.pyplot as plt
import numpy as np
import os


def plot_chemical_abundance(input_dir, output_dir, run, maximum_snap):

    for desired_snap in range(maximum_snap):
        input_path = os.path.join(input_dir, f"snap_{desired_snap:04d}_gas.npz")
        output_path = os.path.join(output_dir, f"MgFe_vs_KFe_{desired_snap:04d}.png")

        data = np.load(input_path)

        ns_Mg = data["ns_Mg"]
        ns_K  = data["ns_K"]
        ns_Fe = data["ns_Fe"]
        os_Mg = data["os_Mg"]
        os_K  = data["os_K"]
        os_Fe = data["os_Fe"]

        ns_x = data["ns_x"]; ns_y = data["ns_y"]; ns_z = data["ns_z"]
        os_x = data["os_x"]; os_y = data["os_y"]; os_z = data["os_z"]

        x  = data["x"];  y  = data["y"];  z  = data["z"]
        mj = data["mj"]

        ns_m = data["ns_m"]
        os_m = data["os_m"]

        total_m = np.sum(mj) + np.sum(ns_m) + np.sum(os_m)

        com_x = (np.sum(x*mj) + np.sum(ns_m*ns_x) + np.sum(os_m*os_x)) / total_m
        com_y = (np.sum(y*mj) + np.sum(ns_m*ns_y) + np.sum(os_m*os_y)) / total_m
        com_z = (np.sum(z*mj) + np.sum(ns_m*ns_z) + np.sum(os_m*os_z)) / total_m

        # keep only stars within 50 pc of COM
        # if 1 sim unit = 100 pc, then 50 pc = 0.5 sim units
        Rcut = 50.0 / 100.0

        r_ns = np.sqrt((ns_x - com_x)**2 + (ns_y - com_y)**2 + (ns_z - com_z)**2)
        r_os = np.sqrt((os_x - com_x)**2 + (os_y - com_y)**2 + (os_z - com_z)**2)

        mask_ns = (r_ns <= Rcut)
        mask_os = (r_os <= Rcut)

        Mg_Fe_ns = (np.log10(ns_Mg/ns_Fe) - np.log10(5.15e-4/1.17e-3))[mask_ns]
        K_Fe_ns  = (np.log10(ns_K/ns_Fe)  - np.log10(3.649e-06/1.17e-3))[mask_ns]

        Mg_Fe_os = (np.log10(os_Mg/os_Fe) - np.log10(5.15e-4/1.17e-3))[mask_os]
        K_Fe_os  = (np.log10(os_K/os_Fe)  - np.log10(3.649e-06/1.17e-3))[mask_os]

        Mg_Fe_obs = (0.20, 0.43, 0.37, 0.42, -0.28, -0.24)
        K_Fe_obs  = (0.29, 0.15, 0.32, 0.27, 0.54, 0.46)

        Mg_Fe_error = [0.04, 0.04, 0.06, 0.05, 0.06, 0.05]
        K_Fe_error  = [0.11, 0.11, 0.12, 0.11, 0.10, 0.12]

        plt.close("all")
        plt.figure()
        plt.scatter(Mg_Fe_ns, K_Fe_ns, s=2, c="red", label="new stars")
        plt.scatter(Mg_Fe_os, K_Fe_os, s=2, c="blue", label="old stars")
        plt.errorbar(Mg_Fe_obs, K_Fe_obs, xerr=Mg_Fe_error, yerr=K_Fe_error, fmt="s", label="observation")
        plt.xlabel("[Mg/Fe]")
        plt.ylabel("[K/Fe]")
        plt.ylim(-0.2, 1.0)
        plt.xlim(-0.5, 0.7)
        plt.tight_layout()
        plt.legend()
        plt.savefig(output_path, dpi=200)
        plt.close()
