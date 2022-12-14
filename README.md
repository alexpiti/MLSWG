# MLSWG
Frequency-domain vector mode-solver for guided waves in **multi-layer slab waveguides (MLSWG)**. 

Can solve for an arbitrary number of slabs (layers), sandwiched between two semi-infinite claddding/substrate layers. Each slab is described by its thickness and refractive index (RI); the RI can be a complex number, allowing the study of lossy and/or plasmonic (SPP) modes. These waveguide structures support only TE- and TM-polarized modes, i.e., neither TEM nor fully hybrid ones. 

![Schematic](https://user-images.githubusercontent.com/97299585/199936285-9da2dab0-859b-4b35-9d20-dac152f73d10.JPG)

<sub>**Fig. 1:** Schematic of a multi-layered slab waveguide (MLSWG), here with two intermediate layers (slabs) between a substrate and cladding (semi-infinite). Propagation is along the z-axis and the 1D cross-section of the waveguide is along the x-axis. The supported mode polarizations are marked in the right.</sub>

## Brief Description
The software works by solving the **characteristic equation** (CE) for a specified wavelength and polarization (TE or TM), and tries to find *all modes* (CE roots) within a specified interval of effective index values; the modal effective index is given by n_eff = β/k0, where β is the propagation constant and k0 is the vacuum wavenumber, both in rad/m. The CE of the MLSWG is formed by applying the EM-field boundary conditions at the interface between each pair of layers; for materials with losses the CE is complex valued. 

The solver returns the eigenvalues (n_effs) and eigenvectors [Ey(x) or Hy(x) profiles, for TE and TM polarizations, respectively] of all modes found. A `leap-frog' algorithm is used to smartly sweep the n_eff search-range provided by the user, looking for roots. A Newton-Raphson method is used to solve the CE in each sub-domain of the search-range. 

## Utilization and I/O arguments

The packages includes the main m-file, **MLSWG.m**, together with some auxiliary m-files:
* ```MLSWG_CharEq.m``` : forms the complex-valued CE (characteristic equation) as a function of n_eff, also complex. So, the roots are found by requiring for CE(n_eff)=0. This function also computes the derivative of the CE w.r.t. the n_eff.
* ```myNewtonRaphson.m``` : [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method) to find **one** x-root of a **complex**-valued function y=f(x), whose x-derivative, df/dx, can also be computed. This methods looks for a root near a supplied (starting guess) point x0, in a x-bounded space.
* ```interpinv.m``` : Graphically find **all** the roots of the **real**-valued function f(x).

To find the waveguide modes, you need to call ```MSLWG.m```, and here is the syntax used:
```
[ neffs , Fy , x_out ] = MLSWG( ModePol , wl , nLR , ns , ts , x , neSL , DoPlotCharEq , DoPlotModeProfiles , DoVocalize )
```
The I/O argument descriptions are as follows:
```
 ==== Outputs ====
  - neffs : vector containing the effective refr indices of modes found
  - Fy    : corresponding mode profiles (Ey or Hy, for TE and TM modes)
  - x_out : cross-section space, used to plot Fy (if input x is not given)

 ==== Inputs / Obligatory ====
  - ModePol : 'TE' or 'TM' (string) -- Polarization of modes
  - wl  : wavelength [Same units as ts!]
  - nLR : refr. indices of [Left,Right] semi-inf layers (2x1 array) -- It's the n_bot and n_top in Fig. 1.
  - ns  : refr. indices of guiding layers (NLx1 array) -- It's the n1 and n2 in Fig. 1.
  - ts  : thicknesses of guiding layers (NLx1 array) [Same units as wl!] -- It's the t1 and t2 in Fig. 1.

 ==== Inputs / Optional ====
  - neSL: n-effective search limits (2x1 array, [Low,High], real valued)
  - x   : WG 1D cross-section (Nx-length vector, real-valued) [Same units as wl & ts!]

 ==== Inputs / Monitoring (optional) ====
  - DoPlotCharEq       : 1 or 0 --> Plot Characteristic Equation
  - DoPlotModeProfiles : 1 or 0 --> Plot mode field profiles found
  - DoVocalize         : 1 or 0 --> Disp info on MATLAB's command window
```
Notes: 
* The cartesian coordinate system (xyz) used follows the convention of **Fig. 1**: The wave propagation direction is along z-axis, the slabs lie in the yz-plane, and the 1D cross-section direction is the x-axis. 
* The first guiding-layer (slab) is assumed to cover x=[0,ts(1)]; the second slab covers x=ts(1)+[0,ts(2)] and so on.

## Examples

**(A) Detuned Photonic Coupler:** In Fig. 2 are the output figures when calling MLSWG with no inputs from the MATLAB command window. In this case, the m-file is set for a lossless photonic coupler with core RI of 3.2 and a substrate/cladding RI of 1.45; the two cores are 250 nm wide and have a gap of 400 nm between them; the operation wavelength is 1550 nm. One of the coupler cores is slightly detuned, i.e., has an RI of 3.18, so that only two quasi-symmetric/antisymmetric modes are supported in the TE polarization.

![OutputFigs_1Default_DesyncCoupler](https://user-images.githubusercontent.com/97299585/199957919-a193aafd-d9af-4a3d-9d2d-9f271cfd678f.JPG)

<sub>**Fig. 2:** Left panels are the CE and its derivative (required for the Newton-Raphson method) in the defined n_eff search-range; the vertical lines correspond to the roots found. The right panels hold the mode profiles with the n_eff given in the panel title.</sub>

**(B) Photonic Lattice:** Below, Fig. 3, an example of the TE modes supported by a lattice of 9 lossy photonic waveguides (same parameters as above, without detuning in any core, but adding a pinch of losses: n=3.2-0.01j).

![OutputFigs_2_Lattice](https://user-images.githubusercontent.com/97299585/199957948-0b86bfff-bbb8-4fb3-a416-ee80981e2f81.JPG)

<sub>**Fig. 3:** The mode profiles supported by a 9-core photonic waveguide lattice.</sub>

**(C) Plasmonic slot:** Finally, in Fig. 4, a gap plasmon mode ([surface plasmon polariton, SPP](https://en.wikipedia.org/wiki/Surface_plasmon_polariton)) supported by a λ-sized air-slot between two semi-infinite gold (Au) layers. The [Au RI in the NIR](https://refractiveindex.info/?shelf=main&book=Au&page=Johnson) is assumed n~0.55-10j, at λ=1550 nm wavelength.

![OutputFigs_3_Plasmonic_GapSPP](https://user-images.githubusercontent.com/97299585/199958064-248dea48-7f4b-4bc3-893c-e20c65227cee.JPG)

<sub>**Fig. 4:** Gap plasmon mode profile (TM-polarization). Plasmonic guiding allows for the reduction of the air-slot width to way below the wavelength, surpassing the diffraction-limit of classical index-constrast waveguiding.</sub>

## Various Useful Notes

* Isotropic materials are mainly considered here, but **anisotropic (birefringent) uniaxial** materials can also be studied provided that their optical axis is aligned with one of the principal cartesian axes (e.g., xyz) and the user carefully defines the RI (ordinary or extraordinary) with respect to the specific mode-polarization (TE or TM). More details can be found [here](https://en.wikipedia.org/wiki/Birefringence).

* The **effective index method** (EIM) can be used, in conjunction with MLSWG, to approximately study the modes in 2D waveguides. More details can be found [here](https://www.computational-photonics.eu/eims.html) (click on `Details' in the top-bar).

* The `leap-frog' algorithm used by MLSWG to **find all roots/modes** is not fully deterministic. So, in some difficult cases (e.g. multi-mode waveguides with multiple and/or electrically thick slabs), the solver might miss some modes; re-executing the code and/or using different neSL inputs helps.

* **Plasmonic (SPP) modes** are trickier for MLSWG to automatically handle, with default inputs only. Do note that TM polarization is specifically required, and a narrow and careful neSL choice is advised. As a last-ditch, edit the MLSWG.m file and switch variable `whichNeffsToOutput' from 0 to 2 (line ~103).

* This solver was partially inspired by Dr. Hammer's [OMS](https://www.computational-photonics.eu/oms.html) which, however, does not support complex-valued refractive indices (nor scripting or editing, of course).

## Citation

If MLSWG was used in research papers, be a pal and cite [this paper of mine](https://doi.org/10.1364/JOSAB.470129) :innocent: EIM and 1D waveguides were preliminarily used there, before switching to full 2D waveguides, solved with the finite-element method.
