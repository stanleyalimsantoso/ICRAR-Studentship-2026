import matplotlib.pyplot as plt
import numpy as np
import os

run = 5

base_dir = r"D:\ICRAR Studentship Data"
base_dir = os.path.join(base_dir, f"Run {run}")

input_path = os.path.join(base_dir, "raw", "sf.dat")
output_path = os.path.join(base_dir, "plots")

with open(input_path, "r") as file:
    lines = file.read().splitlines()
    last_line = lines[-1]
    data = lines[:-1]

    nstep = int(last_line)

    time_list = []
    mass_list = []

    for block in range(nstep):

        time = data[block*3]
        n1, n2, n3 = data[block*3 + 1].split()
        m1, m2, m3 = data[block*3 + 2].split()

        time_list.append(float(time))
        mass_list.append(float(m1))

    time_list = np.array(time_list)
    mass_list = np.array(mass_list)

    time_diff = time_list - np.roll(time_list,1)
    time_diff = time_diff[1:]
    time_diff *= 1.4e8

    mass_diff = mass_list - np.roll(mass_list,1)
    mass_diff = mass_list[1:]
    mass_diff *= 6e10

    SFR_list = mass_diff/time_diff

    plt.plot(time_diff, SFR_list)
    plt.title("Star Formation Rate vs Time")
    plt.xlabel("Time(yr)")
    plt.ylabel("Star Formation Rate")
    output_path = os.path.join(output_path, "SFR vs time.png")
    plt.savefig(output_path)

