# Input audit — 0.9.0 (2026-10-12)

Official references:

- [Poisson solvers](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/POISSON.html)
- [Hybrid k-point calculations](https://manual.cp2k.org/cp2k-2025_1-branch/methods/dft/hartree-fock/ri_kpoints.html)
- [Hybrid interaction potential](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/XC/HF/INTERACTION_POTENTIAL.html)
- [Dispersion parameters](https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/XC/VDW_POTENTIAL/PAIR_POTENTIAL.html)

Partial periodicities now write matching electrostatic boundary conditions:
MT for NONE/XY, ANALYTIC for the other 1D/2D choices. ANALYTIC has slow
convergence, so vacuum/grid convergence remains a calculation requirement.
Hybrid cutoff radii round down to stay inside the inscribed cell sphere.
Numeric settings retain precision. D2 does not request the D3 parameter file.

The writer rejects direct HFX with a k-point mesh; CP2K supports such calculations
through RI-HFX, which this generator does not configure. Smeared diagonalisation
with zero added orbitals is rejected. Preview and save show these errors without
writing invalid input, and recover when corrected. Imports in the plugin entry
point, dialog and writer are at module scope.

Tests cover section structure, all periodicities, numerical precision, cutoff
safety, and Qt preview/save recovery. No CP2K executable is installed here;
documentation and regression checks do not replace an engine calculation.
