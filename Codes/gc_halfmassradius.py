import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree



run = 81
base_dir = rf"D:\ICRAR Studentship Data\Run {run}\GC Data"
preprocessed_dir = os.path.join(base_dir, "preprocessed")
radial_dir = os.path.join(base_dir, "radial")
os.makedirs(radial_dir, exist_ok=True)

max_timestamp = 10

pc_per_unit = 100.0
Rmax_pc = 50.0
Rmax_sim = Rmax_pc / pc_per_unit  # 0.5 sim units

nbins = 500
edges_pc = np.linspace(0.0, Rmax_pc, nbins + 1)
centers_pc = 0.5 * (edges_pc[:-1] + edges_pc[1:] )
shell_vol = (4.0 / 3.0) * np.pi * (edges_pc[1:]**3 - edges_pc[:-1]**3)  # pc^3


def knn_peak_center(px, py, pz, pm, k_nn=64, n_refine=2):
    """
    Pick the particle in the densest region using kNN:
    densest = smallest distance to k-th nearest neighbor.

    Returns (cx, cy, cz).
    """
    pos = np.column_stack([px, py, pz])
    n = pos.shape[0]
    if n == 0:
        return np.nan, np.nan, np.nan

    kk = min(k_nn + 1, n)  # +1 because the nearest neighbor is itself
    tree = cKDTree(pos)
    dists, _ = tree.query(pos, k=kk)
    rk = dists[:, -1]

    i0 = int(np.argmin(rk))
    c = pos[i0].copy()

    r_ref = rk[i0]
    if (not np.isfinite(r_ref)) or r_ref <= 0.0:
        return c[0], c[1], c[2]

    for _ in range(n_refine):
        ids = tree.query_ball_point(c, r_ref)
        if len(ids) < 1:
            break
        w = pm[ids]
        wsum = np.sum(w)
        if (not np.isfinite(wsum)) or wsum <= 0.0:
            break
        c = np.sum(pos[ids] * w[:, None], axis=0) / wsum

    return c[0], c[1], c[2]


def shrinking_sphere_center(sx, sy, sz, sm, cx0, cy0, cz0, n_iter=30, shrink=0.8, min_keep=50):
    """
    Same shrinking-sphere idea as your current code, but you can pass an initial guess.
    """
    cx, cy, cz = cx0, cy0, cz0

    for _ in range(n_iter):
        dx = sx - cx
        dy = sy - cy
        dz = sz - cz
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        rmax = np.max(r)

        if (not np.isfinite(rmax)) or rmax <= 0.0:
            break

        keep = r < shrink * rmax
        if keep.sum() < min_keep:
            break

        sx2 = sx[keep]; sy2 = sy[keep]; sz2 = sz[keep]; sm2 = sm[keep]
        mtot2 = np.sum(sm2)
        if (not np.isfinite(mtot2)) or mtot2 <= 0.0:
            break

        cx = np.sum(sm2 * sx2) / mtot2
        cy = np.sum(sm2 * sy2) / mtot2
        cz = np.sum(sm2 * sz2) / mtot2

    return cx, cy, cz


# Tunables
k_nn = 64
Rlocal_pc = 50.0
Rlocal_sim = Rlocal_pc / pc_per_unit

for desired_snap in range(0, max_timestamp):
    input_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")
    data = np.load(input_path)

    # gas (optional, not used for centering here)
    x = data["x"]; y = data["y"]; z = data["z"]; mj = data["mj"]

    # stars
    ns_x = data["ns_x"]; ns_y = data["ns_y"]; ns_z = data["ns_z"]; ns_m = data["ns_m"]
    os_x = data["os_x"]; os_y = data["os_y"]; os_z = data["os_z"]; os_m = data["os_m"]

    # build star arrays for centering (sim units, mass units)
    sx = np.concatenate([ns_x, os_x])
    sy = np.concatenate([ns_y, os_y])
    sz = np.concatenate([ns_z, os_z])
    sm = np.concatenate([ns_m, os_m])

    if sm.size == 0 or np.sum(sm) <= 0.0:
        print(f"snap {desired_snap:04d}: no stars at all")
        continue

    # --- NEW: density-peak centre via KD-tree kNN, then local shrinking-sphere ---
    # Use new stars for peak-finding if there are enough of them, otherwise fall back to all stars.
    if ns_x.size >= (k_nn + 1):
        cx0, cy0, cz0 = knn_peak_center(ns_x, ns_y, ns_z, ns_m, k_nn=k_nn, n_refine=2)
    else:
        cx0, cy0, cz0 = knn_peak_center(sx, sy, sz, sm, k_nn=min(k_nn, max(1, sx.size - 1)), n_refine=2)

    # Restrict to a local region around the peak before shrinking-sphere
    dx0 = sx - cx0; dy0 = sy - cy0; dz0 = sz - cz0
    r0 = np.sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0)
    local = r0 < Rlocal_sim

    if local.sum() >= 50:
        sxL = sx[local]; syL = sy[local]; szL = sz[local]; smL = sm[local]
        cx, cy, cz = shrinking_sphere_center(sxL, syL, szL, smL, cx0, cy0, cz0, n_iter=30, shrink=0.8, min_keep=50)
    else:
        # Not enough stars in the local cut, just use the kNN peak centre
        cx, cy, cz = cx0, cy0, cz0

    # recenter everything on density centre
    x = x - cx; y = y - cy; z = z - cz
    ns_x = ns_x - cx; ns_y = ns_y - cy; ns_z = ns_z - cz
    os_x = os_x - cx; os_y = os_y - cy; os_z = os_z - cz

    # radii in pc
    r_ns_pc = np.sqrt(ns_x*ns_x + ns_y*ns_y + ns_z*ns_z) * pc_per_unit
    r_os_pc = np.sqrt(os_x*os_x + os_y*os_y + os_z*os_z) * pc_per_unit

    # restrict to within 50 pc for Kenji’s request
    ns_in = r_ns_pc < Rmax_pc
    os_in = r_os_pc < Rmax_pc

    r_ns_pc_50 = r_ns_pc[ns_in]
    r_os_pc_50 = r_os_pc[os_in]
    m_ns_msun_50 = ns_m[ns_in] * 1e6
    m_os_msun_50 = os_m[os_in] * 1e6

    r_star_pc_50 = np.concatenate([r_ns_pc_50, r_os_pc_50])
    m_star_msun_50 = np.concatenate([m_ns_msun_50, m_os_msun_50])

    if r_star_pc_50.size == 0:
        print(f"snap {desired_snap:04d}: no stars within 50 pc after centering")
        continue

    # M(<R): sort by radius then cumulative sum
    order = np.argsort(r_star_pc_50)
    r_sorted = r_star_pc_50[order]
    m_sorted = m_star_msun_50[order]
    Menc = np.cumsum(m_sorted)

    Mtot50 = Menc[-1]
    target = 0.5 * Mtot50

    k = int(np.searchsorted(Menc, target, side="left"))
    if k == 0:
        R_half_pc = r_sorted[0]
    else:
        M1 = Menc[k - 1]; M2 = Menc[k]
        R1 = r_sorted[k - 1]; R2 = r_sorted[k]
        if M2 == M1:
            R_half_pc = R2
        else:
            R_half_pc = R1 + (target - M1) / (M2 - M1) * (R2 - R1)

    print(
        f"snap {desired_snap:04d}: "
        f"M(<50 pc) = {Mtot50:.6e} Msun, "
        f"half-mass radius = {R_half_pc:.6f} pc, "
        f"Nstars(<50pc) = {r_star_pc_50.size}"
    )

    plt.figure()
    if r_sorted.size == 1:
        plt.scatter(r_sorted, Menc)
    else:
        plt.plot(r_sorted, Menc)
    plt.xlim(0.0, Rmax_pc)
    plt.xlabel("R (pc)")
    plt.ylabel("M(<R) (Msun)")
    plt.savefig(
        os.path.join(radial_dir, f"Mcum_stars_within50pc_snap_{desired_snap:04d}.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # radial density profile rho(r) in spherical shells (within 50 pc)
    hist_ns, _ = np.histogram(r_ns_pc_50, bins=edges_pc, weights=m_ns_msun_50)
    hist_os, _ = np.histogram(r_os_pc_50, bins=edges_pc, weights=m_os_msun_50)
    hist_tot = hist_ns + hist_os

    rho_ns = np.full(nbins, np.nan)
    rho_os = np.full(nbins, np.nan)
    rho_tot = np.full(nbins, np.nan)

    np.divide(hist_ns, shell_vol, out=rho_ns, where=(shell_vol > 0.0))
    np.divide(hist_os, shell_vol, out=rho_os, where=(shell_vol > 0.0))
    np.divide(hist_tot, shell_vol, out=rho_tot, where=(shell_vol > 0.0))

    plt.figure()
    plt.plot(centers_pc, rho_tot, label="Total stars")
    plt.plot(centers_pc, rho_ns, label="New stars")
    plt.plot(centers_pc, rho_os, label="Old stars")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(0.0, Rmax_pc)
    plt.xlabel("R (pc)")
    plt.ylabel("rho (Msun/pc^3)")
    plt.legend()
    plt.savefig(
        os.path.join(radial_dir, f"rho_stars_within20pc_snap_{desired_snap:04d}.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
