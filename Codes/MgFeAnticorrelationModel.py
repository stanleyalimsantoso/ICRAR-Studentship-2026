import numpy as np
import matplotlib.pyplot as plt

mg_Fe_obs = [0.20, 0.43, 0.37, 0.42, -0.28, -0.24]
K_Fe_obs = [0.29, 0.15, 0.32, 0.27, 0.54, 0.46]

change_y = False

MGMC = 1

sAGB = 0 #0 for illiadis (K/Fe = 1.5 and Mg/Fe = -1.5),
         #1 for ventura (K/Fe = 0.6, Mg/Fe = -0.2)

R_K_GMC = 10**(0) #This is the ratio of the chemical abundance of [K/Fe] in GMCs
R_Mg_GMC = 10**(0.5) #This is the ratio of the chemical abundance of [Mg/Fe] in GMCs

R_K_Massive = 10**(0) #This is the ratio of the chemical abundance of [K/Fe] in Massive AGB stars (5-7 solar mass)
R_Mg_Massive = 10**(-0.2) #This is the ratio of the chemical abundance of [Mg/Fe] in Massive AGB stars (5-7 solar mass)

K_mix_list = []
Mg_mix_list = []

if sAGB == 0: 
    R_K_Super = 10**(1.5) #This is the ratio of the chemical abundance of [K/Fe] in Super AGB stars (7-10 solar mass)
    R_Mg_Super = 10**(-1.5) #This is the ratio of the chemical abundance of [Mg/Fe] in Super AGB stars (7-10 solar mass)
    tag = "illiadis"

elif sAGB ==1: 
    R_K_Super = 10**0.6
    R_Mg_Super = 10**(-0.2)
    tag = "ventura"


if change_y == True:
    Msuper = np.linspace(0,100,600)
    Mmassive = 50

    for mej in Msuper:
        R_K_mix = (MGMC*R_K_GMC + mej*R_K_Super + Mmassive*R_K_Massive)/(MGMC+mej+Mmassive)
        R_Mg_mix = (MGMC*R_Mg_GMC + mej*R_Mg_Super + Mmassive*R_Mg_Massive)/(MGMC+mej+Mmassive)

        K_mix_list.append(np.log10(R_K_mix))
        Mg_mix_list.append(np.log10(R_Mg_mix))

    plt.scatter(K_mix_list, Mg_mix_list, c="blue", s=5)
    plt.plot(K_mix_list, Mg_mix_list)
    plt.scatter(K_Fe_obs, mg_Fe_obs, c="red", s=5)
    plt.ylabel("[Mg/Fe]")
    plt.xlabel("[K/Fe]")
    plt.title("[K/Fe] vs [Mg/Fe]")
    plt.savefig(f"D:/ICRAR Studentship Data/KMg plot models/z = {Mmassive}, {tag}.jpg")
    plt.show()

else:
    Mmassive = np.linspace(0,5,50)
    Msuper = 10

    for mej in Mmassive:
        R_K_mix = (MGMC*R_K_GMC + Msuper*R_K_Super + mej*R_K_Massive)/(MGMC+Msuper+mej)
        R_Mg_mix = (MGMC*R_Mg_GMC + Msuper*R_Mg_Super + mej*R_Mg_Massive)/(MGMC+Msuper+mej)

        K_mix_list.append(np.log10(R_K_mix))
        Mg_mix_list.append(np.log10(R_Mg_mix))

    plt.scatter(K_mix_list, Mg_mix_list, c="blue", s=5)
    plt.plot(K_mix_list, Mg_mix_list)
    plt.scatter(K_Fe_obs, mg_Fe_obs, c="red", s=5)
    plt.ylabel("[Mg/Fe]")
    plt.xlabel("[K/Fe]")
    plt.title("[K/Fe] vs [Mg/Fe]")
    plt.savefig(f"D:/ICRAR Studentship Data/KMg plot models/z = {Msuper}, {tag}.jpg")
    plt.show()

