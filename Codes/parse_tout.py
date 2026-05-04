import os
import numpy as np
import argparse

from classify_gmc import classify_gmc
from newstar_dist import newstar_dist
from data_output import data_output
from heatmap import heatmap

#rsync -avh --progress salimsan@ozstar.swin.edu.au:/fred/oz080/runh2m1.dir/tout.dat "/mnt/d/ICRAR Studentship Data/Run 16/"
#rsync -avh --progress salimsan@ozstar.swin.edu.au:/fred/oz080/runvimf300.dir/tout.dat "/mnt/d/ICRAR Studentship Data/Run 11/GC data"
# rsync -avPxH --no-g --chmod=Dg+s "/mnt/d/ICRAR Studentship Data/Run 60/GC Data/raw/tout.dat" salimsan@data-mover01.hpc.swin.edu.au:/home/salimsan/runvimf300.dir/gmc.dat
def parse_tout(input_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    current_snap = 0

    with open(input_path, "r") as file:
        while True:
            header = file.readline().split()
            if not header:
                break

            num_bodies = int(header[0])
            current_time=float(header[1])

            max_snap = int(header[2])

            file.readline()
            file.readline()
            file.readline()

            x, y, z, vx, vy, vz, mj, frah, t = [], [], [], [], [], [], [], [], []
            ns_x, ns_y, ns_z, ns_t, ns_m, ns_vx, ns_vy, ns_vz = [], [], [], [], [], [], [], []
            os_m, dm_m = [], []

            for i in range(num_bodies):
                line = file.readline().split()
                if len(line) < 14:
                    continue

                iwas = int(float(line[6]))
                if iwas == 2:
                    x.append(float(line[0]))
                    y.append(float(line[1]))
                    z.append(float(line[2]))
                    vx.append(float(line[3]))
                    vy.append(float(line[4]))
                    vz.append(float(line[5]))
                    mj.append(float(line[8]))
                    t.append(float(line[12]))
                    frah.append(float(line[13]))

                elif iwas == 3:
                    ns_x.append(float(line[0]))
                    ns_y.append(float(line[1]))
                    ns_z.append(float(line[2]))
                    ns_vx.append(float(line[3]))
                    ns_vy.append(float(line[4]))
                    ns_vz.append(float(line[5]))
                    ns_t.append(float(line[12]))
                    ns_m.append(float(line[8]))

                elif iwas == 1:
                    os_m.append(float(line[8]))
                elif iwas ==0:
                    dm_m.append(float(line[8]))

            x = np.array(x)
            y = np.array(y)
            z = np.array(z)
            vx = np.array(vx)
            vy = np.array(vy)
            vz = np.array(vz)
            mj = np.array(mj)
            frah = np.array(frah)
            t = np.array(t)

            ns_x = np.array(ns_x)
            ns_y = np.array(ns_y)
            ns_z = np.array(ns_z)
            ns_vx = np.array(ns_vx)
            ns_vy = np.array(ns_vy)
            ns_vz = np.array(ns_vz) 
            ns_t = np.array(ns_t)
            ns_m = np.array(ns_m)

            os_m = np.array(os_m)
            dm_m = np.array(dm_m)


            np.savez_compressed(
                f"{output_dir}/snap_{current_snap:04d}_gas.npz",
                x=x, y=y, z=z,
                vx=vx, vy=vy, vz=vz,
                mj=mj, frah=frah, t=t,
                ns_x=ns_x, ns_y=ns_y, ns_z=ns_z, ns_t = ns_t, ns_m = ns_m,
                ns_vx = ns_vx, ns_vy = ns_vy, ns_vz = ns_vz, os_m = os_m, dm_m = dm_m
            )

            print("Saved snapshot", current_snap, "with", len(x), "gas particles")

            if current_snap > max_snap:
                break

            current_snap += 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--preprocessed", required=True)
    p.add_argument("--max-timestamp", type=int, required=True)
    p.add_argument("--total-gas", action="store_true")
    p.add_argument("--gmc-dir", required=True)
    p.add_argument("--run", type=int, required=True)
    p.add_argument("--postprocessed", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--plots-dir", required=True)



    args = p.parse_args()
    parse_tout(args.input, args.preprocessed)
    classify_gmc(args.preprocessed, args.gmc_dir, args.max_timestamp, args.total_gas, args.run)
    newstar_dist(args.gmc_dir, args.postprocessed, args.max_timestamp, args.run, args.total_gas)
    data_output(args.postprocessed, args.output, args.max_timestamp, args.run, args.total_gas )
    heatmap(args.gmc_dir, args.plots_dir, args.run, args.max_timestamp, args.total_gas)

if __name__ == "__main__":
    main()