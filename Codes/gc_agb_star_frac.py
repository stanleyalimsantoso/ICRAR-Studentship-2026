import os
import numpy as np

max_timestamp = 10
run = 59

base_path = f"D:\\ICRAR Studentship Data\\Run {run}\\GC Data"
preprocessed_dir = os.path.join(base_path, "preprocessed")

t_thr_upper = 4/3
t_thr_lower = 2/3
r_thr = 2

for desired_snap in range(max_timestamp):
    input_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")
    data = np.load(input_path)

    # Load data
    ns_t = data["ns_t"]
    ns_m = data["ns_m"]
    ns_x = data["ns_x"]
    ns_y = data["ns_y"]
    ns_z = data["ns_z"]

    os_m = data["os_m"]
    os_x = data["os_x"]
    os_y = data["os_y"]
    os_z = data["os_z"]

    # Calculate center of mass for old stars
    os_total_mass = os_m.sum()
    os_com_x = np.sum(os_m * os_x) / os_total_mass
    os_com_y = np.sum(os_m * os_y) / os_total_mass
    os_com_z = np.sum(os_m * os_z) / os_total_mass

    # Calculate center of mass for new stars
    ns_total_mass = ns_m.sum()
    ns_com_x = np.sum(ns_m * ns_x) / ns_total_mass
    ns_com_y = np.sum(ns_m * ns_y) / ns_total_mass
    ns_com_z = np.sum(ns_m * ns_z) / ns_total_mass

    # Combine to get overall center of mass
    total_mass = os_total_mass + ns_total_mass
    com_x = (np.sum(os_m * os_x) + np.sum(ns_m * ns_x)) / total_mass
    com_y = (np.sum(os_m * os_y) + np.sum(ns_m * ns_y)) / total_mass
    com_z = (np.sum(os_m * os_z) + np.sum(ns_m * ns_z)) / total_mass

    # Renormalize coordinates based on overall COM
    os_x_com = os_x - com_x
    os_y_com = os_y - com_y
    os_z_com = os_z - com_z
    
    ns_x_com = ns_x - com_x
    ns_y_com = ns_y - com_y
    ns_z_com = ns_z - com_z

    # Now use these renormalized coordinates for your analysis
    # Create masks with renormalized coordinates

    os_r_thr = np.sqrt(os_x_com**2 + os_y_com**2 + os_z_com**2)
    ns_r_thr = np.sqrt(ns_x_com**2 + ns_y_com**2 + ns_z_com**2)
    
    agb_mask = (ns_t < t_thr_upper) & (ns_t > t_thr_lower)
    os_mask = (os_r_thr < r_thr)
    ns_mask = (ns_r_thr < r_thr)


    # Apply masks correctly - combine AGB and GC masks

    
    # Calculate masses
    total_gc_mass = os_m[os_mask].sum() + ns_m[ns_mask].sum()
    agb_mass = ns_m[agb_mask].sum()
    
    # Avoid division by zero
    if total_gc_mass > 0:
        frac = agb_mass / total_gc_mass
    else:
        frac = 0
    
    print(f"  AGB fraction in GC (r < 0.5 from COM) = {frac:.4f}, timestep = {desired_snap}")
    print("-" * 50)
