# MLSWG
Frequency-domain vector mode solver for guided waves in slab (1D) waveguides. Supports arbitrary number of slabs (layers) sandwiched between two semi-infinite claddding/substrate layers. Each slab is described by its thickness and refractive index; the latter can be a complex number, allowing the study of lossy and/or plasmonic (SPP) modes. These waveguide structures support only TE- and TM-polarized modes, i.e., neither TEM nor fully hybrid ones. 

## Brief Description
The software works by solving the characteristic equation (CE) for the specified wavelength and polarization, and tries to find all modes (CE roots) within a given interval of effective index values (n_eff = β/k0). The CE is formed by applying the EM-field boundary conditions at the interface between each pair of layers; for materials with losses the CE is complex valued. 

The solver returns the eigenvalues (n_effs) and eigenvectors (Ey or Hy profile, for TE and TM polarizations, respectively) of all modes found. A `leap-frog' algorithm is used to smartly sweep the n_eff search-range boundary provided by the user, looking for roots. A Newton-Raphson method is used to solve the CE in each sub-domain of the search-range. 

## Utilization and I/O arguments

The packages includes the main m-file (MLSWG) together with some more auxiliary m-files:
* MLSWG_CharEq.m : forms the CE (characteristic equation).
* myNewtonRaphson.m : [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method) to find **one** x-root of a complex-valued function y=f(x), given also its x-derivative, df/dx. This methods looks for a root near a supplied (starting guess) point x0, in a x-bounded space.
* interpinv.m : Graphically find **all** the roots of the real-valued function f(x).

To find your modes, you need to call MSLWG. Here is the syntax used:

[ neffs , Fy , x_out ] = MLSWG( ModePol , wl , nLR , ns , ts , x , neSL , DoPlotCharEq , DoPlotModeProfiles , DoVocalize )

Here are the I/O argument descriptions:

 ==== Outputs ====
  - neffs : vector containing the effective refr indices of modes found
  - Fy    : corresponding mode profiles (Ey or Hy, for TE and TM modes)
  - x_out : cross-section space, used to plot Fy (if input x is not given)

 ==== Inputs / Obligatory ====
  - ModePol : 'TE' or 'TM' (string) - Polarization of modes
  - wl  : wavelength [Same units as ts!]
  - nLR : refr. indices of [Left,Right] semi-inf layers (2x1 array)
  - ns  : refr. indices of guiding layers (NLx1 array)
  - ts  : thicknesses of guiding layers (NLx1 array) [Same units as wl!]

 ==== Inputs / Optional ====
  - neSL: n-effective search limits (2x1 array, [Low,High])
  - x   : WG 1D cross-section (Nx-length vector) [Same units as wl & ts!]

 ==== Inputs / Monitoring (optional) ====
  - DoPlotCharEq       : 1 or 0 --> Plot Characteristic Equation
  - DoPlotModeProfiles : 1 or 0 --> Plot mode field profiles found
  - DoVocalize         : 1 or 0 --> Disp info on MATLAB's command window

Note: The wave propagation direction is along z-axis and the slabs lie in the yz-plane, so that the 1D cross-section direction is the x-axis. The first guiding-layer (slab) is assumed to cover x=[0,ts(1)]; the second slab covers x=ts(1)+[0,ts(2)] and so on.

## Examples

In Fig. 1 are the output figures when calling MLSWG with no inputs from the MATLAB command window. In this case, the m-file is set for a lossless photonic coupler with core index 3.2 and substrate index 1.45; the two cores are 250 nm wide and have a gap of 400 nm between them; the operation wavelength is 1550 nm. One of the coupler cores is slightly detuned, i.e., has an index of 3.1750, so that only two quasi-symmetric/antisymmetric modes are supported in the TE polarization.

![OutputFigs_1Default_DesyncCoupler](https://user-images.githubusercontent.com/97299585/199725650-5983ebb8-ba11-4e1c-a531-251bc15c552d.JPG)

Fig. 1: Left panels are the CE and its derivative (required for the Newton-Raphson method) in the defined n_eff search-range; the vertical lines correspond to the roots found. The right panels hold the mode profiles with the n_eff given in the panel title.

Below, Fig. 2, an example of the modes supported by a lattice of 9 photonic waveguides (parameters same as above, without detuning in any core).

![OutputFigs_2_Lattice](https://user-images.githubusercontent.com/97299585/199725661-82325518-7c27-4b81-b5c1-1392d9a7b8c5.JPG)

Fig. 2: The mode profiles supported by a 9-core photonic waveguide lattice.

Finally, a gap plasmon mode, when a 250 nm air-slot is formed between two gold layers (Au refractive index in the NIR is assumed n~0.55-10j for 1550 nm wavelength).

![OutputFigs_3_Plasmonic_GapSPP](https://user-images.githubusercontent.com/97299585/199728604-21b66072-732c-45f3-9c4a-a4193ebba6de.JPG)

Fig. 3: Gap plasmon mode profile (TM-polarization).

## Various Useful Notes

* Isotropic materials are only considered, but anisotropic uniaxial ones can also be studied provided that their optical axis is aligned with one of the principal cartesian ones and the user carefully specifies the refractive index (ordinary or extraordinary) for the specific mode-polarization (TE or TM).

* The effective index method (EIM) can be used, in conjunction with MLSWG, to approximately study the modes in 2D waveguides. More details can be found [here](https://www.computational-photonics.eu/eims.html) (click on `Details' in the top-bar).

* As the `leap-frog' algorithm used to find all roots/modes is not fully deterministic, in some difficult cases (e.g. multi-mode waveguides with multiple and/or thick slabs), the solver might miss some modes; rerunning the code and/or using different neSL inputs helps.

* Plasmonic (SPP) modes are trickier for MLSWG to automatically handle. Do note that TM polarization is specificially required, and a narrow and careful neSL choice is advised. As a last-ditch, switch whichNeffsToOutput from 0 to 2 (line ~97 in the MLSWG m-file).

* Partially inspired by Dr. Hammer's [OMS](https://www.computational-photonics.eu/oms.html), which does not support complex-valued refractive indices (nor scripting or editing, of course).
