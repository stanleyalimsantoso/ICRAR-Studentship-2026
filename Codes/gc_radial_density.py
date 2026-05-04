import os
import numpy as np
import matplotlib.pyplot as plt

max_timestamp = 5
run = 30

base_path = f"D:\ICRAR Studentship Data\Run {run}\GC Data"
preprocessed_dir = os.path.join(base_path, "preprocessed")
radial_dir = os.path.join(base_path, "radial")

os.makedirs(radial_dir, exist_ok=True)

bins = 300

for desired_snap in range(max_timestamp):
    snap_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")
    data = np.load(snap_path)

    x = data["x"]; y = data["y"]; z = data["z"]; mj = data["mj"]
    ns_x = data["ns_x"]; ns_y = data["ns_y"]; ns_z = data["ns_z"]; ns_m = data["ns_m"]
    os_x = data["os_x"]; os_y = data["os_y"]; os_z = data["os_z"]; os_m = data["os_m"]

    total_m = np.sum(mj) + np.sum(ns_m) + np.sum(os_m)

    com_x = (np.sum(x*mj) + np.sum(ns_m*ns_x) + np.sum(os_m*os_x)) / total_m
    com_y = (np.sum(y*mj) + np.sum(ns_m*ns_y) + np.sum(os_m*os_y)) / total_m
    com_z = (np.sum(z*mj) + np.sum(ns_m*ns_z) + np.sum(os_m*os_z)) / total_m

    x = x - com_x; y = y - com_y; z = z - com_z
    ns_x = ns_x - com_x; ns_y = ns_y - com_y; ns_z = ns_z - com_z
    os_x = os_x - com_x; os_y = os_y - com_y; os_z = os_z - com_z

    r_thr_lower = 0.5
    r_thr_upper = 1.0

    os_r2 = os_x**2 + os_y**2 + os_z**2
    ns_r2 = ns_x**2 + ns_y**2 + ns_z**2

    oslowermask = (os_r2 < r_thr_lower**2)
    osuppermask = (os_r2 < r_thr_upper**2)

    nslowermask = (ns_r2 < r_thr_lower**2)
    nsuppermask = (ns_r2 < r_thr_upper**2)

    m_os_lower = os_m[oslowermask]
    m_os_upper = os_m[osuppermask]
    m_ns_lower = ns_m[nslowermask]
    m_ns_upper = ns_m[nsuppermask]

    out_txt = os.path.join(radial_dir, f"mass_in_spheres_snap_{desired_snap:04d}.txt")

    vals = np.array([
        m_os_lower.sum()*10**6,
        m_os_upper.sum()*10**6,
        m_ns_lower.sum()*10**6,
        m_ns_upper.sum()*10**6,
    ], dtype=float)

    np.savetxt(
        out_txt,
        vals.reshape(1, -1),   # 1 row, 4 columns
        fmt="%.8e",
        header="OS_lower_mass_sum  OS_upper_mass_sum  NS_lower_mass_sum  NS_upper_mass_sum"
    )

    for tag, xx, yy, mm in [
        ("normal", x, y, mj),
        ("ns", ns_x, ns_y, ns_m),
        ("os", os_x, os_y, os_m),
    ]:
        plt.figure()
        plt.gcf().set_facecolor("black")
        plt.gca().set_facecolor("black")

        if xx.size > 0:
            Lx = np.percentile(np.abs(xx), 99.5)
            Ly = np.percentile(np.abs(yy), 99.5)
            L = max(Lx, Ly)
        else:
            L = 1.0

        if (not np.isfinite(L)) or (L <= 0.0):
            L = 1.0

        H, xedges, yedges = np.histogram2d(
            xx, yy,
            bins=bins,
            range=[[-L, L], [-L, L]],
            weights=mm
        )

        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        Sigma = H / (dx * dy)

        logSigma = np.full_like(Sigma, np.nan, dtype=float)
        mask = Sigma > 0.0
        if np.any(mask):
            logSigma[mask] = np.log10(Sigma[mask])
            lo = np.nanpercentile(logSigma[mask], 5)
            hi = np.nanpercentile(logSigma[mask], 99.5)
            if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                img = (np.clip(logSigma, lo, hi) - lo) / (hi - lo)
            else:
                img = np.zeros_like(Sigma, dtype=float)
        else:
            img = np.zeros_like(Sigma, dtype=float)

        im = plt.imshow(
            img.T,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            cmap="plasma",
            aspect="equal",
            interpolation="nearest",
            vmin=0.0,
            vmax=1.0
        )
        plt.colorbar(im)

        out_png = os.path.join(radial_dir, f"{tag}_snap_{desired_snap:04d}.png")
        plt.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close()

    Rmax = 1
    nbins = 20
    edges = np.linspace(0.0, Rmax, nbins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    r_ns = np.sqrt(ns_x*ns_x + ns_y*ns_y + ns_z*ns_z)
    r_os = np.sqrt(os_x*os_x + os_y*os_y + os_z*os_z)

    n_ns_bin, _ = np.histogram(r_ns, bins=edges)
    n_os_bin, _ = np.histogram(r_os, bins=edges)

    denom = n_ns_bin + n_os_bin
    frac_new = np.zeros_like(denom, dtype=float)
    np.divide(n_ns_bin, denom, out=frac_new, where=(denom > 0))

    plt.figure()
    plt.plot(centers, frac_new)
    plt.ylim(0.0, 1.0)
    plt.xlabel("R (sim units)")   # 0..0.00571 corresponds to 0..100 pc
    plt.ylabel("New-star fraction (count)")
    out_png = os.path.join(radial_dir, f"frac_new_count_snap_{desired_snap:04d}.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    os_Fe = data["os_Fe"]
    os_K = data["os_K"]
    os_Mg = data["os_Mg"]

    ns_Fe = data["ns_Fe"]
    ns_K = data["ns_K"]
    ns_Mg = data["ns_Mg"]

    Rmax_sim = 1.0
    Rmax_pc = 100.0
    nbins = 20

    edges_sim = np.linspace(0.0, Rmax_sim, nbins + 1)
    centers_sim = 0.5 * (edges_sim[:-1] + edges_sim[1:])
    centers_pc = centers_sim * Rmax_pc

    r_ns = np.sqrt(ns_x*ns_x + ns_y*ns_y + ns_z*ns_z)
    r_os = np.sqrt(os_x*os_x + os_y*os_y + os_z*os_z)

    Fe_sun = 1.17e-3
    Mg_sun = 5.15e-4
    K_sun  = 3.94e-7

    ns_mask = (ns_Fe > 0) & (ns_K > 0) & (ns_Mg > 0)
    os_mask = (os_Fe > 0) & (os_K > 0) & (os_Mg > 0)

    kfe_ns = np.full(ns_Fe.shape, np.nan)
    mgfe_ns = np.full(ns_Fe.shape, np.nan)
    kfe_os = np.full(os_Fe.shape, np.nan)
    mgfe_os = np.full(os_Fe.shape, np.nan)

    kfe_ns[ns_mask]  = np.log10(ns_K[ns_mask]  / ns_Fe[ns_mask]) - np.log10(K_sun  / Fe_sun)
    mgfe_ns[ns_mask] = np.log10(ns_Mg[ns_mask] / ns_Fe[ns_mask]) - np.log10(Mg_sun / Fe_sun)

    kfe_os[os_mask]  = np.log10(os_K[os_mask]  / os_Fe[os_mask]) - np.log10(K_sun  / Fe_sun)
    mgfe_os[os_mask] = np.log10(os_Mg[os_mask] / os_Fe[os_mask]) - np.log10(Mg_sun / Fe_sun)

    r_all = np.concatenate([r_ns[ns_mask], r_os[os_mask]])
    kfe_all = np.concatenate([kfe_ns[ns_mask], kfe_os[os_mask]])
    mgfe_all = np.concatenate([mgfe_ns[ns_mask], mgfe_os[os_mask]])

    m_all = np.concatenate([ns_m[ns_mask], os_m[os_mask]]) * 1e6  # Msun

    m_bin, _   = np.histogram(r_all, bins=edges_sim, weights=m_all)
    mk_bin, _  = np.histogram(r_all, bins=edges_sim, weights=m_all * kfe_all)
    mmg_bin, _ = np.histogram(r_all, bins=edges_sim, weights=m_all * mgfe_all)

    mean_kfe  = np.full(nbins, np.nan)
    mean_mgfe = np.full(nbins, np.nan)
    np.divide(mk_bin,  m_bin, out=mean_kfe,  where=(m_bin > 0.0))
    np.divide(mmg_bin, m_bin, out=mean_mgfe, where=(m_bin > 0.0))

    plt.figure()
    plt.plot(centers_pc, mean_kfe)
    plt.xlabel("R (pc)")
    plt.ylabel("Mean [K/Fe]")
    plt.savefig(os.path.join(radial_dir, f"mean_KFe_snap_{desired_snap:04d}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.plot(centers_pc, mean_mgfe)
    plt.xlabel("R (pc)")
    plt.ylabel("Mean [Mg/Fe]")
    plt.savefig(os.path.join(radial_dir, f"mean_MgFe_snap_{desired_snap:04d}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()


    # --- stellar (new + old) radial mass density profile rho(r) ---
    # convention: 1 sim unit = 100 pc, mass unit = 1e6 Msun

    r_ns = np.sqrt(ns_x*ns_x + ns_y*ns_y + ns_z*ns_z)   # sim units
    r_os = np.sqrt(os_x*os_x + os_y*os_y + os_z*os_z)   # sim units

    m_ns = ns_m * 1e6   # Msun
    m_os = os_m * 1e6   # Msun

    # consistent with centers_pc = centers_sim * 100 (your code)
    edges_pc = edges_sim * Rmax_pc  # pc

    m_ns_bin, _ = np.histogram(r_ns, bins=edges_sim, weights=m_ns)
    m_os_bin, _ = np.histogram(r_os, bins=edges_sim, weights=m_os)

    shell_vol = (4.0/3.0) * np.pi * (edges_pc[1:]**3 - edges_pc[:-1]**3)   # pc^3

    rho_ns = np.full(nbins, np.nan)
    rho_os = np.full(nbins, np.nan)
    rho_star = np.full(nbins, np.nan)

    np.divide(m_ns_bin, shell_vol, out=rho_ns, where=(shell_vol > 0.0))                 # Msun/pc^3
    np.divide(m_os_bin, shell_vol, out=rho_os, where=(shell_vol > 0.0))                 # Msun/pc^3
    np.divide(m_ns_bin + m_os_bin, shell_vol, out=rho_star, where=(shell_vol > 0.0))    # Msun/pc^3

    plt.figure()
    plt.plot(centers_pc, rho_ns, label="New stars")
    plt.plot(centers_pc, rho_os, label="Old stars")
    plt.plot(centers_pc, rho_star, label="Total stars")
    plt.yscale("log")
    plt.xlabel("R (pc)")
    plt.ylabel("Stellar density rho (Msun/pc^3)")
    plt.legend()
    plt.savefig(os.path.join(radial_dir, f"rho_stars_snap_{desired_snap:04d}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
