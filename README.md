# ICRAR-Studentship-2026
This repository contains an overview of my Summer Studentship project at ICRAR, in collaboration with Pawsey. It includes the code and results for the project, titled “Computer Simulations of Globular Cluster Formation in the Early LMC,” supervised by Kenji Bekki.

This repository is still a work in progress. I’ve attached my final report, and you can also see a short video about my project on the Pawsey alumni page (check the 2025–2026 interns): https://pawsey.org.au/alumni/

## What I did

The project is about simulating the formation of globular clusters in the early Large Magellanic Cloud, with a focus on NGC 1786. This cluster is interesting because recent observations suggest that it shows a rare Mg-K anticorrelation, where some stars are magnesium-poor but potassium-rich.

The main question is whether this chemical pattern can be reproduced using computer simulations. The idea I am testing is that massive or super-AGB stars may have polluted nearby gas with material that is low in magnesium and high in potassium. If later stars formed from different mixtures of this polluted gas and normal gas, this could explain the Mg-K anticorrelation.

In this repo, I will collect the code, figures, and results from the project. The goal is not just to make a simulation that looks nice, but to find model globular clusters that have reasonable physical properties, such as mass and half-mass radius, while also reproducing the observed chemical abundance pattern.

<p align="center">
  <img src="pictures/Picture1.png" alt="Mg-K anti-correlation plot" width="550">
</p>

The physical mechanism behind the Mg-K anti-correlation is still not fully understood. One difficulty is that this pattern is rare, and there are only a small number of observed stars and clusters where this effect has been measured clearly. Because of this, simulations are useful because they let us test whether different enrichment scenarios can actually reproduce the observed chemical abundance pattern.

One possible explanation is pollution from massive or super-AGB stars. Super-AGB stars can reach very high internal temperatures, where magnesium can be destroyed and potassium can be produced through nuclear reactions. As a result, the material they eject is expected to be relatively low in Mg and high in K. These stars then lose this processed material through stellar winds, which can mix with the surrounding pristine gas inside or near the forming globular cluster.

If new stars form from gas with different mixtures of pristine material and super-AGB ejecta, they can naturally end up with different Mg and K abundances. Stars forming from more polluted gas would appear Mg-poor and K-rich, while stars forming from less polluted gas would look more chemically normal. This gives a possible way to explain the observed Mg-K anti-correlation in NGC 1786.

<p align="center">
  <img src="pictures/Picture5.png" alt="visualisation of the sAGB polluter model" width="550">
</p>

But first what are super AGB stars? AGB stars are evolved stars that have already used up the hydrogen in their core. At this stage, the star has a compact core, surrounded by helium-burning and hydrogen-burning shells, with a large outer envelope. The core no longer produces energy in the same way as a normal main-sequence star, so most of the energy comes from these burning shells. The outer envelope is cool, extended, and not strongly bound, which means the star can lose a lot of material through stellar winds.

Massive AGB stars are simply the higher-mass end of AGB stars. Because they are more massive, the base of their outer envelope can become hot enough for extra nuclear processing, called hot bottom burning. This is important because the material in the envelope can be chemically changed before it is lost through stellar winds. Super-AGB stars are even more massive, and they can ignite carbon in their core, leaving behind an oxygen-neon core instead of the usual carbon-oxygen core. For this project, these stars matter because their winds can carry chemically processed material into nearby gas, which may later form new stars with different Mg and K abundances.

Super-AGB stars have very hot hydrogen-burning regions, especially through hot bottom burning, where the base of the convective envelope becomes hot enough for proton-capture reactions. In these conditions, magnesium nuclei can capture protons and be converted into heavier elements through the Mg-Al cycle, so the remaining gas becomes lower in magnesium. At even higher temperatures, proton captures involving argon can produce potassium, so the material expelled by the star can become potassium-rich.

<p align="center">
  <img src="pictures/Picture4.png" alt="diagram of the core of AGB stars" width="550">
</p>

Super-AGB stars are useful for this problem because they can reach very high temperatures during hot bottom burning. In this stage, the base of the convective envelope becomes hot enough for proton-capture reactions. Magnesium can be reduced through the Mg-Al cycle, while potassium can be enhanced through high-temperature reactions involving argon. This means the material expelled by super-AGB stars can become Mg-poor and K-rich.

In the polluter model, this processed material is expelled through stellar winds and mixes with the surrounding pristine gas in the forming globular cluster. Stars that form mostly from pristine gas should have more normal Mg and K abundances, while stars that form from gas with more super-AGB ejecta should become more Mg-poor and K-rich. This gives a simple way to produce an Mg-K anti-correlation.

As a first-pass analytical check, I used a two-step polluter model. The first step mixes pristine gas with super-AGB ejecta, which gives the strong Mg-poor and K-rich enrichment needed for the most extreme stars. However, super-AGB pollution alone does not smoothly explain all the observed stars, so the second step adds massive AGB ejecta. Massive AGB ejecta is less extreme, so it can adjust the abundance pattern without pushing the gas too far toward very low Mg and very high K.

The abundance of each element is calculated in linear space first:

$$
A_X = \frac{A_{X,0} + f_{\mathrm{sAGB}}A_{X,\mathrm{sAGB}} + f_{\mathrm{mAGB}}A_{X,\mathrm{mAGB}}}{1 + f_{\mathrm{sAGB}} + f_{\mathrm{mAGB}}}
$$

where

$$
A_X = 10^{[X/\mathrm{Fe}]}
$$

After mixing, the abundance is converted back into bracket notation:

$$
[X/\mathrm{Fe}]_{\mathrm{mix}} = \log_{10}(A_X)
$$

Here, \(f_{\mathrm{sAGB}}\) and \(f_{\mathrm{mAGB}}\) are the mass ratios of super-AGB and massive AGB ejecta compared to the original GMC gas mass.

| Gas source | Role in the model | [Mg/Fe] yield | [K/Fe] yield | Meaning |
|---|---|---:|---:|---|
| Pristine gas | Starting gas | 0.0 | 0.0 | Normal gas before AGB pollution |
| Super-AGB ejecta | Main polluter | -1.5 | 1.5 | Strongly Mg-poor and K-rich material |
| Massive AGB ejecta | Secondary polluter | -0.3 | 0.0 | Mildly Mg-poor material used to adjust the model |

Using these assumed yields, the model gives the following approximate fits:

| ID | sAGB ejecta mass ratio | mAGB ejecta mass ratio | [Mg/Fe] | [K/Fe] |
|---:|---:|---:|---:|---:|
| 978  | 0.080 | 1.520 | 0.20  | 0.29 |
| 1321 | 0.016 | 0.207 | 0.43  | 0.15 |
| 1436 | 0.052 | 0.407 | 0.37  | 0.32 |
| 1501 | 0.035 | 0.220 | 0.92  | 0.27 |
| 2310 | 4.220 | 54.600 | -0.28 | 0.54 |
| 2418 | 4.220 | 54.600 | -0.24 | 0.46 |

The first four stars can be explained reasonably well with the super-AGB pollution model. The last two stars can also be reproduced by the two-step model, but only by using a very high massive-AGB-to-GMC mass ratio. This is probably not physically ideal, but it is still useful because it shows that the polluter model is worth investigating further. It also suggests that the assumed yields may need to be improved or fine-tuned.



