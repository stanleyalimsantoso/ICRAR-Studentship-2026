import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

max_timestamp = 10
run = 81

base_path = f"D:\\ICRAR Studentship Data\\Run {run}\\GC Data"
preprocessed_dir = os.path.join(base_path, "preprocessed")
plots_dir = os.path.join(base_path, "plots")
os.makedirs(plots_dir, exist_ok=True)

R_thr = 5.0  # sim units

SOL_MG_FE = 5.15e-4 / 1.17e-3
SOL_K_FE  = 3.649e-06 / 1.17e-3

# eyeballed chemical boxes (dex)
# (name, Mg_min, Mg_max, K_min, K_max)
boxes = [
    ("A", 0.44, 0.62, -0.12, 0.06),
    ("B", 0.42, 0.52,  0.08, 0.18),
    ("C", 0.32, 0.44,  0.22, 0.33),
    ("D", 0.14, 0.33,  0.33, 0.47),
    ("E",-0.03, 0.05,  0.48, 0.57),
    ("F",-0.42,-0.34,  0.56, 0.64),
    ("G",-0.26,-0.16,  0.56, 0.62),
]

# fixed colors for each island (use the same mapping in every plot)
_TAB = plt.get_cmap("tab10")
GROUP_COLORS = {name: _TAB((i + 1) % 10) for i, (name, _, _, _, _) in enumerate(boxes)}



# surface-density knobs
SIM_TO_PC = 100.0
BIN_PC = 0.5
PAD_PC = 5.0
OVERLAY_ALPHA = 0.65

for desired_snap in range(max_timestamp):

    input_path = os.path.join(preprocessed_dir, f"snap_{desired_snap:04d}_gas.npz")
    data = np.load(input_path)

    ns_x  = data["ns_x"]
    ns_y  = data["ns_y"]
    ns_z  = data["ns_z"]
    ns_m  = data["ns_m"]
    ns_K  = data["ns_K"]
    ns_Mg = data["ns_Mg"]
    ns_Fe = data["ns_Fe"]

    total_m = np.sum(ns_m)
    if total_m <= 0:
        continue

    # COM shift
    com_x = np.sum(ns_x * ns_m) / total_m
    com_y = np.sum(ns_y * ns_m) / total_m
    com_z = np.sum(ns_z * ns_m) / total_m
    ns_x = ns_x - com_x
    ns_y = ns_y - com_y
    ns_z = ns_z - com_z

    # spatial cut
    GC_mask = (ns_x**2 + ns_y**2 + ns_z**2) <= R_thr**2

    x  = ns_x[GC_mask]
    y  = ns_y[GC_mask]
    m  = ns_m[GC_mask]        # <-- this was missing in your file
    Mg = ns_Mg[GC_mask]
    K  = ns_K[GC_mask]
    Fe = ns_Fe[GC_mask]

    # valid chemistry
    valid = (Mg > 0) & (K > 0) & (Fe > 0) & (m > 0)
    x  = x[valid]
    y  = y[valid]
    m  = m[valid]
    Mg = Mg[valid]
    K  = K[valid]
    Fe = Fe[valid]

    if len(m) == 0:
        continue

    Mg_Fe = np.log10(Mg / Fe) - np.log10(SOL_MG_FE)
    K_Fe  = np.log10(K  / Fe) - np.log10(SOL_K_FE)

    # label by boxes
    labels = np.full(len(Mg_Fe), "noise", dtype=object)
    for name, mg0, mg1, k0, k1 in boxes:
        sel = (Mg_Fe >= mg0) & (Mg_Fe <= mg1) & (K_Fe >= k0) & (K_Fe <= k1)
        labels[sel] = name

    # pc grid
    x_pc = x * SIM_TO_PC
    y_pc = y * SIM_TO_PC
    m = m*10**10
    x_edges = np.arange(x_pc.min() - PAD_PC, x_pc.max() + PAD_PC + BIN_PC, BIN_PC)
    y_edges = np.arange(y_pc.min() - PAD_PC, y_pc.max() + PAD_PC + BIN_PC, BIN_PC)

    # background surface density
    H_tot, _, _ = np.histogram2d(x_pc, y_pc, bins=[x_edges, y_edges], weights=m)
    Sigma_tot = H_tot / (BIN_PC * BIN_PC)
    logSigma_tot = np.log10(Sigma_tot + 1e-30)

    vmin = np.percentile(logSigma_tot, 5)
    vmax = np.percentile(logSigma_tot, 99)
        # --- focus window: 150 x 150 pc square centered on (0,0) ---
    FOCUS_PC = 300.0
    HALF = FOCUS_PC / 2.0

    # crop particles to the window
    in_win = (np.abs(x_pc) <= HALF) & (np.abs(y_pc) <= HALF)
    x_pc = x_pc[in_win]
    y_pc = y_pc[in_win]
    m    = m[in_win]
    labels = labels[in_win]

    # fixed histogram edges for the same 150x150 square
    x_edges = np.arange(-HALF, HALF + BIN_PC, BIN_PC)
    y_edges = np.arange(-HALF, HALF + BIN_PC, BIN_PC)

    extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]


    plt.figure(figsize=(8, 8))
    plt.imshow(np.zeros_like(logSigma_tot.T), origin="lower", extent=extent, cmap="gray", vmin=0, vmax=1)


    # colored overlays by island
    handles = []

    for i, (name, _, _, _, _) in enumerate(boxes):
        sel = (labels == name)
        if not np.any(sel):
            continue

        H_r, _, _ = np.histogram2d(x_pc[sel], y_pc[sel], bins=[x_edges, y_edges], weights=m[sel])
        Sigma_r = H_r / (BIN_PC * BIN_PC)
        log_r = np.log10(Sigma_r + 1e-30)

        a = (log_r - vmin) / (vmax - vmin)
        a = np.clip(a, 0.0, 1.0)
        a[Sigma_r <= 0] = 0.0
        a = OVERLAY_ALPHA * a

        col = GROUP_COLORS[name]

        rgba = np.zeros((a.T.shape[0], a.T.shape[1], 4))
        rgba[..., 0] = col[0]
        rgba[..., 1] = col[1]
        rgba[..., 2] = col[2]
        rgba[..., 3] = a.T

        plt.imshow(rgba, origin="lower", extent=extent, interpolation="nearest")
        handles.append(Patch(facecolor=col, edgecolor="none", label=name))

    plt.gca().set_aspect("equal", adjustable="box")
    from matplotlib.ticker import FuncFormatter

    ax = plt.gca()
    x0_pc = com_x * SIM_TO_PC
    y0_pc = com_y * SIM_TO_PC

    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v + x0_pc:.0f}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f"{v + y0_pc:.0f}"))

    ax.set_xlabel("x [pc]")
    ax.set_ylabel("y [pc]")

    if handles:
        plt.legend(handles=handles, frameon=True, loc="upper right")

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, f"xy_surface_islands_{desired_snap:04d}.png"), dpi=200)
    plt.close()
