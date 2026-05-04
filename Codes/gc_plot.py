import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.patches import Circle

try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None


def knn_peak_center(px, py, pz, pm, k_nn=64, n_refine=2):
    """
    Pick the densest particle using kNN: smallest distance to k-th nearest neighbor.
    Then do a couple of mass-weighted refinements inside that local scale.
    """
    if cKDTree is None:
        raise ImportError("scipy is required for cKDTree. Install with: pip install scipy")

    pos = np.column_stack([px, py, pz])
    n = pos.shape[0]
    if n == 0:
        return np.nan, np.nan, np.nan

    kk = min(k_nn + 1, n)  # +1 because nearest is itself
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
    Shrinking-sphere refinement starting from (cx0, cy0, cz0).
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

        sx2 = sx[keep]
        sy2 = sy[keep]
        sz2 = sz[keep]
        sm2 = sm[keep]

        mtot2 = np.sum(sm2)
        if (not np.isfinite(mtot2)) or mtot2 <= 0.0:
            break

        cx = np.sum(sm2 * sx2) / mtot2
        cy = np.sum(sm2 * sy2) / mtot2
        cz = np.sum(sm2 * sz2) / mtot2

    return cx, cy, cz


def gc_plot(input_dir, plots_dir, run, max_timestamp):

    # unit convention: 1 simulation unit = 100 pc
    pc_per_unit = 100.0

    # radii are now all in pc
    Rmax_pc = 20.0
    R_circle = Rmax_pc  # in pc

    # centering knobs (pc)
    k_nn = 64
    Rlocal_pc = 20.0

    for desired_snap in range(max_timestamp):

        gc_path = os.path.join(input_dir, f"snap_{desired_snap:04d}_gas.npz")

        if not os.path.exists(gc_path):
            print(f"Snapshot {desired_snap}: {gc_path} not found, skipping.")
            continue

        data = np.load(gc_path)

        # new stars for the image (convert to pc)
        ns_x = data["ns_x"] * pc_per_unit
        ns_y = data["ns_y"] * pc_per_unit
        ns_m = data["ns_m"] * 10**6
        ns_z = data["ns_z"] * pc_per_unit

        # old stars for centering (if present) (convert to pc)
        if "os_x" in data.files and "os_y" in data.files and "os_m" in data.files:
            os_x = data["os_x"] * pc_per_unit
            os_y = data["os_y"] * pc_per_unit
            os_m = data["os_m"] * 10**6
            os_z = data["os_z"] * pc_per_unit
        else:
            os_x = np.array([], dtype=float)
            os_y = np.array([], dtype=float)
            os_z = np.array([], dtype=float)
            os_m = np.array([], dtype=float)

        # build arrays for centering (now in pc)
        sx = np.concatenate([ns_x, os_x])
        sy = np.concatenate([ns_y, os_y])
        sz = np.concatenate([ns_z, os_z])
        sm = np.concatenate([ns_m, os_m])

        if sm.size == 0 or np.sum(sm) <= 0.0:
            cx, cy, cz = np.nan, np.nan, np.nan
        else:
            # kNN peak centre: prefer ns if enough points, else all stars
            if ns_x.size >= (k_nn + 1):
                cx0, cy0, cz0 = knn_peak_center(ns_x, ns_y, ns_z, ns_m, k_nn=k_nn, n_refine=2)
            else:
                kk = min(k_nn, max(1, sx.size - 1))
                cx0, cy0, cz0 = knn_peak_center(sx, sy, sz, sm, k_nn=kk, n_refine=2)

            # local cut then shrink refine (pc)
            dx0 = sx - cx0
            dy0 = sy - cy0
            dz0 = sz - cz0
            r0 = np.sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0)
            local = r0 < Rlocal_pc

            if local.sum() >= 50:
                cx, cy, cz = shrinking_sphere_center(
                    sx[local], sy[local], sz[local], sm[local],
                    cx0, cy0, cz0,
                    n_iter=30, shrink=0.8, min_keep=50
                )
            else:
                cx, cy, cz = cx0, cy0, cz0

        # choose your plotting radius (pc)
        Rplot_pc = 20.0

        if np.isfinite(cx) and np.isfinite(cy) and np.isfinite(cz):
            dxp = ns_x - cx
            dyp = ns_y - cy
            dzp = ns_z - cz
            keep_plot = (dxp*dxp + dyp*dyp + dzp*dzp) < (Rplot_pc * Rplot_pc)

            px = ns_x[keep_plot]
            py = ns_y[keep_plot]
            pm = ns_m[keep_plot]
        else:
            px, py, pm = ns_x, ns_y, ns_m

        # ---- plotting (same structure as your original) ----
        plt.figure()
        plt.gcf().set_facecolor("black")
        ax = plt.gca()
        ax.set_facecolor("black")

        bins = 300

        H, xedges, yedges = np.histogram2d(
            px, py,
            bins=bins,
            range=[[cx - Rplot_pc, cx + Rplot_pc],
                   [cy - Rplot_pc, cy + Rplot_pc]],
            weights=pm
        )

        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        Sigma = H / (dx * dy)  # now mass / pc^2

        logSigma = np.full_like(Sigma, np.nan, dtype=float)
        mask = Sigma > 0.0
        if np.any(mask):
            logSigma[mask] = np.log10(Sigma[mask])
            lo = np.nanpercentile(logSigma[mask], 0)
            hi = np.nanpercentile(logSigma[mask], 100)

            if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
                img = np.clip(logSigma, lo, hi)
                vmin, vmax = lo, hi
            else:
                vmin = np.nanmin(logSigma[mask])
                vmax = np.nanmax(logSigma[mask])
                img = logSigma
        else:
            img = np.zeros_like(Sigma, dtype=float)
            vmin, vmax = 0.0, 1.0

        im = plt.imshow(
            img.T,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            cmap="plasma",
            aspect="equal",
            interpolation="nearest",
            vmin=vmin,
            vmax=vmax
        )

        cbar = plt.colorbar(im)
        cbar.set_label(r"$\log_{10}(\Sigma)\ [\mathrm{Msun}\,\mathrm{pc}^{-2}]$", color="white")
        cbar.ax.yaxis.set_tick_params(color="white")
        plt.setp(plt.getp(cbar.ax.axes, "yticklabels"), color="white")

        # optional axis labels (now meaningful in pc)
        ax.set_xlabel("x [pc]", color="white")
        ax.set_ylabel("y [pc]", color="white")
        ax.tick_params(colors="white")

        os.makedirs(plots_dir, exist_ok=True)

        out_png = os.path.join(plots_dir, f"snap_{desired_snap:04d}.png")
        plt.savefig(out_png, dpi=300, bbox_inches="tight")

        # overlay centre and circle (pc), then save second figure
        if np.isfinite(cx) and np.isfinite(cy):
            ax.plot([cx], [cy], marker="x", markersize=7, markeredgewidth=1.2,
                    color="white", linestyle="None")
            circ = Circle((cx, cy), R_circle, fill=False, edgecolor="white",
                          linewidth=1.0, alpha=0.9)
            ax.add_patch(circ)

        out_png2 = os.path.join(plots_dir, f"snap_{desired_snap:04d}_center50pc.png")
        plt.savefig(out_png2, dpi=300, bbox_inches="tight")

        plt.close()
