import os
import numpy as np
import argparse


from gc_plot import gc_plot
from plot_chemical_abundance import plot_chemical_abundance

def parse_tout_gc(input_path, output_dir, max_timestamp):

    os.makedirs(output_dir, exist_ok=True)

    with open(input_path, "r") as file:
        current_snap = 0
        header = file.readline().split()
        num_bodies = int(float(header[0]))
        placeholder_line = file.readline().split()

        while current_snap < max_timestamp:

            while True:
                raw = file.readline()
                if raw == "":
                    return  # EOF
                line = raw.split()
                if not line:
                    continue
                if len(line) < 5:
                    num_bodies = int(float(line[0]))
                    break

            x, y, z, vx, vy, vz, mj, frah, t = [], [], [], [], [], [], [], [], []
            ns_x, ns_y, ns_z, ns_t, ns_m, ns_vx, ns_vy, ns_vz = [], [], [], [], [], [], [], []
            os_x, os_y, os_z, os_t, os_m, os_vx, os_vy, os_vz = [], [], [], [], [], [], [], []

            ns_H, ns_He, ns_C, ns_N, ns_O, ns_Fe, ns_Mg, ns_Ca, ns_Si, ns_Na, ns_Ba, ns_Eu, ns_Al, ns_K = [], [], [], [], [], [], [], [], [], [], [], [], [], []
            os_H, os_He, os_C, os_N, os_O, os_Fe, os_Mg, os_Ca, os_Si, os_Na, os_Ba, os_Eu, os_Al, os_K = [], [], [], [], [], [], [], [], [], [], [], [], [], []
            gas_H, gas_He, gas_C, gas_N, gas_O, gas_Fe, gas_Mg, gas_Ca, gas_Si, gas_Na, gas_Ba, gas_Eu, gas_Al, gas_K = [], [], [], [], [], [], [], [], [], [], [], [], [], []     
            
            for i in range(num_bodies):
                line = file.readline().split()
                if not line:
                    break

                if len(line) < 30:
                    continue

                iwas = int(float(line[6]))

                if iwas == 2:
                    x.append(float(line[0]))
                    y.append(float(line[1]))
                    z.append(float(line[2]))
                    vx.append(float(line[3]))
                    vy.append(float(line[4]))
                    vz.append(float(line[5]))

                    mj.append(float(line[10]))
                    t.append(float(line[14]))
                    frah.append(float(line[15]))

                    gas_H.append(float(line[16]))
                    gas_He.append(float(line[17]))
                    gas_C.append(float(line[18]))
                    gas_N.append(float(line[19]))
                    gas_O.append(float(line[20]))
                    gas_Fe.append(float(line[21]))
                    gas_Mg.append(float(line[22]))
                    gas_Ca.append(float(line[23]))
                    gas_Si.append(float(line[24]))
                    gas_Na.append(float(line[25]))
                    gas_Ba.append(float(line[26]))
                    gas_Eu.append(float(line[27]))
                    gas_Al.append(float(line[28]))
                    gas_K.append(float(line[29]))

                elif iwas == 3:
                    ns_x.append(float(line[0]))
                    ns_y.append(float(line[1]))
                    ns_z.append(float(line[2]))
                    ns_vx.append(float(line[3]))
                    ns_vy.append(float(line[4]))
                    ns_vz.append(float(line[5]))

                    ns_t.append(float(line[14]))    # tt (formation time)
                    ns_m.append(float(line[10]))    # mj

                    ns_H.append(float(line[16]))
                    ns_He.append(float(line[17]))
                    ns_C.append(float(line[18]))
                    ns_N.append(float(line[19]))
                    ns_O.append(float(line[20]))
                    ns_Fe.append(float(line[21]))
                    ns_Mg.append(float(line[22]))
                    ns_Ca.append(float(line[23]))
                    ns_Si.append(float(line[24]))
                    ns_Na.append(float(line[25]))
                    ns_Ba.append(float(line[26]))
                    ns_Eu.append(float(line[27]))
                    ns_Al.append(float(line[28]))
                    ns_K.append(float(line[29]))

                elif iwas == 1:
                    os_x.append(float(line[0]))
                    os_y.append(float(line[1]))
                    os_z.append(float(line[2]))
                    os_vx.append(float(line[3]))
                    os_vy.append(float(line[4]))
                    os_vz.append(float(line[5]))

                    os_t.append(float(line[14]))    # tt (formation time)
                    os_m.append(float(line[10]))    # mj

                    os_H.append(float(line[16]))
                    os_He.append(float(line[17]))
                    os_C.append(float(line[18]))
                    os_N.append(float(line[19]))
                    os_O.append(float(line[20]))
                    os_Fe.append(float(line[21]))
                    os_Mg.append(float(line[22]))
                    os_Ca.append(float(line[23]))
                    os_Si.append(float(line[24]))
                    os_Na.append(float(line[25]))
                    os_Ba.append(float(line[26]))
                    os_Eu.append(float(line[27]))
                    os_Al.append(float(line[28]))
                    os_K.append(float(line[29]))

            x = np.array(x); y = np.array(y); z = np.array(z)
            vx = np.array(vx); vy = np.array(vy); vz = np.array(vz)
            mj = np.array(mj); frah = np.array(frah); t = np.array(t)

            gas_H = np.array(gas_H); gas_He = np.array(gas_He); gas_C = np.array(gas_C); gas_N = np.array(gas_N)
            gas_O = np.array(gas_O); gas_Fe = np.array(gas_Fe); gas_Mg = np.array(gas_Mg); gas_Ca = np.array(gas_Ca)
            gas_Si = np.array(gas_Si); gas_Na = np.array(gas_Na); gas_Ba = np.array(gas_Ba); gas_Eu = np.array(gas_Eu)
            gas_Al = np.array(gas_Al); gas_K = np.array(gas_K)

            ns_x = np.array(ns_x); ns_y = np.array(ns_y); ns_z = np.array(ns_z)
            ns_vx = np.array(ns_vx); ns_vy = np.array(ns_vy); ns_vz = np.array(ns_vz)
            ns_t = np.array(ns_t); ns_m = np.array(ns_m)

            ns_H = np.array(ns_H); ns_He = np.array(ns_He); ns_C = np.array(ns_C); ns_N = np.array(ns_N)
            ns_O = np.array(ns_O); ns_Fe = np.array(ns_Fe); ns_Mg = np.array(ns_Mg); ns_Ca = np.array(ns_Ca)
            ns_Si = np.array(ns_Si); ns_Na = np.array(ns_Na); ns_Ba = np.array(ns_Ba); ns_Eu = np.array(ns_Eu)
            ns_Al = np.array(ns_Al); ns_K = np.array(ns_K)

            os_x = np.array(os_x); os_y = np.array(os_y); os_z = np.array(os_z)
            os_vx = np.array(os_vx); os_vy = np.array(os_vy); os_vz = np.array(os_vz)
            os_t = np.array(os_t); os_m = np.array(os_m)

            os_H = np.array(os_H); os_He = np.array(os_He); os_C = np.array(os_C); os_N = np.array(os_N)
            os_O = np.array(os_O); os_Fe = np.array(os_Fe); os_Mg = np.array(os_Mg); os_Ca = np.array(os_Ca)
            os_Si = np.array(os_Si); os_Na = np.array(os_Na); os_Ba = np.array(os_Ba); os_Eu = np.array(os_Eu)
            os_Al = np.array(os_Al); os_K = np.array(os_K)

            np.savez_compressed(
                f"{output_dir}/snap_{current_snap:04d}_gas.npz",
                x=x, y=y, z=z,
                vx=vx, vy=vy, vz=vz,
                mj=mj, frah=frah, t=t,
                gas_H=gas_H, gas_He=gas_He, gas_C=gas_C, gas_N=gas_N,
                gas_O=gas_O, gas_Fe=gas_Fe, gas_Mg=gas_Mg, gas_Ca=gas_Ca,
                gas_Si=gas_Si, gas_Na=gas_Na, gas_Ba=gas_Ba, gas_Eu=gas_Eu,
                gas_Al=gas_Al, gas_K=gas_K,
                ns_x=ns_x, ns_y=ns_y, ns_z=ns_z, ns_t=ns_t, ns_m=ns_m,
                ns_vx=ns_vx, ns_vy=ns_vy, ns_vz=ns_vz,
                ns_H=ns_H, ns_He=ns_He, ns_C=ns_C, ns_N=ns_N,
                ns_O=ns_O, ns_Fe=ns_Fe, ns_Mg=ns_Mg, ns_Ca=ns_Ca,
                ns_Si=ns_Si, ns_Na=ns_Na, ns_Ba=ns_Ba, ns_Eu=ns_Eu,
                ns_Al=ns_Al, ns_K=ns_K,
                os_x=os_x, os_y=os_y, os_z=os_z, os_t=os_t, os_m=os_m,
                os_vx=os_vx, os_vy=os_vy, os_vz=os_vz,
                os_H=os_H, os_He=os_He, os_C=os_C, os_N=os_N,
                os_O=os_O, os_Fe=os_Fe, os_Mg=os_Mg, os_Ca=os_Ca,
                os_Si=os_Si, os_Na=os_Na, os_Ba=os_Ba, os_Eu=os_Eu,
                os_Al=os_Al, os_K=os_K
            )

            print("Saved snapshot", current_snap, "with", len(x), "gas particles")
            current_snap += 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--preprocessed", required=True)
    p.add_argument("--max-timestamp", type=int, required=True)
    p.add_argument("--run", type=int, required=True)
    p.add_argument("--plots-dir", required=True)


    args = p.parse_args()
    parse_tout_gc(args.input, args.preprocessed, args.max_timestamp)
    plot_chemical_abundance(args.preprocessed, args.plots_dir, args.run, args.max_timestamp)
    gc_plot(args.preprocessed, args.plots_dir, args.run, args.max_timestamp)


if __name__ == "__main__":
    main()