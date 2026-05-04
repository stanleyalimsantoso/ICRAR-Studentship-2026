import numpy as np
from scipy.spatial import cKDTree

def naive_search(x, y, z, vx, vy, vz, mj, frah, r_thr, m_thr, total_gas):
    if total_gas == True:
        frah = np.ones(len(frah))
    
    thr2 = r_thr * r_thr
    ngas = len(x)
    visited = bytearray(ngas)

    gmc_x_com = []
    gmc_y_com = []
    gmc_z_com = []
    gmc_m_tot = []
    gmc_n_part = []

    labels = -1 * np.ones(ngas, dtype=int)
    gmc_id = 0

    for i in range(ngas):
        if visited[i]:
            continue

        members = [i]
        visited[i] = 1
        k = 0

        while k < len(members):
            a = members[k]
            k += 1
            xa, ya, za = x[a], y[a], z[a]

            for j in range(ngas):
                if visited[j]:
                    continue
                dx = xa - x[j]
                dy = ya - y[j]
                dz = za - z[j]
                if dx*dx + dy*dy + dz*dz <= thr2:
                    visited[j] = 1
                    members.append(j)

        if len(members) < 2:
            continue

        m_tot = 0.0
        x_tot = y_tot = z_tot = 0.0

        for j in members:
            m = mj[j] * frah[j]
            m_tot += m
            x_tot += m * x[j]
            y_tot += m * y[j]
            z_tot += m * z[j]

        if m_tot < m_thr:
            continue

        gmc_x_com.append(x_tot / m_tot)
        gmc_y_com.append(y_tot / m_tot)
        gmc_z_com.append(z_tot / m_tot)
        gmc_m_tot.append(m_tot)
        gmc_n_part.append(len(members))

        for j in members:
            labels[j] = gmc_id

        gmc_id += 1

    gmc_x_com = np.array(gmc_x_com)
    gmc_y_com = np.array(gmc_y_com)
    gmc_z_com = np.array(gmc_z_com)
    gmc_m_tot = np.array(gmc_m_tot)
    gmc_n_part = np.array(gmc_n_part)

    return gmc_x_com, gmc_y_com, gmc_z_com, gmc_m_tot, gmc_n_part, labels


def cdktree(x, y, z, vx, vy, vz, mj, frah, r_thr, m_thr_upper, m_thr_lower, total_gas):
    linking_length = 0.0005
    if total_gas == True:
        frah = np.ones(len(frah))

    pos = np.vstack((x, y, z)).T
    ngas = len(x)
    tree = cKDTree(pos)
    visited = np.zeros(ngas, dtype=bool)

    gmc_x_com = []
    gmc_y_com = []
    gmc_z_com = []
    gmc_vx_com = []
    gmc_vy_com = []
    gmc_vz_com = []
    gmc_m_tot = []
    gmc_n_part = []
    potential_gmc_mass = []

    labels = -1 * np.ones(ngas, dtype=int)
    gmc_id = 0

    for i in range(ngas):
        if visited[i]:
            continue

        members = [i]
        visited[i] = True
        k = 0

        while k < len(members):
            a = members[k]
            k += 1
            neighbours = tree.query_ball_point(pos[a], linking_length)
            for j in neighbours:
                if not visited[j]:
                    visited[j] = True
                    members.append(j)

        if len(members) < 2:
            continue

        m_tot = 0.0
        x_tot = y_tot = z_tot = 0.0
        vx_tot = vy_tot = vz_tot = 0.0

        for j in members:
            m = mj[j] * frah[j]
            m_tot += m
            x_tot += m * x[j]
            y_tot += m * y[j]
            z_tot += m * z[j]
            vx_tot += m * vx[j]
            vy_tot += m * vy[j]
            vz_tot += m * vz[j]

        if m_tot == 0:
            continue

        x_com = x_tot / m_tot
        y_com = y_tot / m_tot
        z_com = z_tot / m_tot

        r2_thr = r_thr * r_thr

        m_tot = 0.0
        x_tot = y_tot = z_tot = 0.0
        vx_tot = vy_tot = vz_tot = 0.0
        kept = []

        for j in members:
            dx = x[j] - x_com
            dy = y[j] - y_com
            dz = z[j] - z_com
            if dx*dx + dy*dy + dz*dz <= r2_thr:
                kept.append(j)
                m = mj[j] * frah[j]
                m_tot += m
                x_tot += m * x[j]
                y_tot += m * y[j]
                z_tot += m * z[j]
                vx_tot += m * vx[j]
                vy_tot += m * vy[j]
                vz_tot += m * vz[j]

        if len(kept) < 2:
            continue

        if m_tot > m_thr_upper:
            potential_gmc_mass.append(m_tot)
            continue
        elif m_tot < m_thr_lower:
            potential_gmc_mass.append(m_tot)

        gmc_x_com.append(x_tot / m_tot)
        gmc_y_com.append(y_tot / m_tot)
        gmc_z_com.append(z_tot / m_tot)
        gmc_vx_com.append(vx_tot / m_tot)
        gmc_vy_com.append(vy_tot / m_tot)
        gmc_vz_com.append(vz_tot / m_tot)
        gmc_m_tot.append(m_tot)
        gmc_n_part.append(len(kept))

        for j in kept:
            labels[j] = gmc_id

        gmc_id += 1

    gmc_x_com = np.array(gmc_x_com)
    gmc_y_com = np.array(gmc_y_com)
    gmc_z_com = np.array(gmc_z_com)
    gmc_vx_com = np.array(gmc_vx_com)
    gmc_vy_com = np.array(gmc_vy_com)
    gmc_vz_com = np.array(gmc_vz_com)
    gmc_m_tot = np.array(gmc_m_tot)
    gmc_n_part = np.array(gmc_n_part)
    potential_gmc_mass = np.array(potential_gmc_mass)

    if len(gmc_m_tot) != 0:
        print(f"Max GMC mass: {np.max(gmc_m_tot)}")

    if len(potential_gmc_mass) != 0:
        print(np.max(potential_gmc_mass))

    pos_lbl = labels[labels >= 0]
    if pos_lbl.size > 0:
        assert pos_lbl.max() < len(gmc_x_com), (pos_lbl.min(), pos_lbl.max(), len(gmc_x_com))

    return gmc_x_com, gmc_y_com, gmc_z_com, gmc_vx_com, gmc_vy_com, gmc_vz_com, gmc_m_tot, gmc_n_part, labels
