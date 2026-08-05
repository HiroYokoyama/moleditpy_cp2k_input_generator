# MoleditPy CP2K Input Generator

A [MoleditPy](https://github.com/HiroYokoyama/python_molecular_editor) plugin that
writes **CP2K (Quickstep)** input files from a molecule or a crystal structure,
with a live preview of the complete input before you save.

## What it does

- **Three structure sources**
  - the molecule currently open in MoleditPy, wrapped in an orthorhombic vacuum box
  - a `.cif` file loaded from disk (asymmetric unit expanded with the CIF's own
    symmetry operations)
  - the structure already open in the **CIF Viewer** plugin panel, copied across
    without re-reading the file
- **Supercells** — independent a/b/c repeats applied to any source
- **Surface slab builder** — cut any (hkl) surface from a bulk structure: layers,
  vacuum thickness, a termination shift that slides the cut window through the
  cell, and an option to put **c** along the surface normal. Hexagonal
  Miller-Bravais indices `(hkil)` are accepted and folded to `(hkl)`
- **Vacuum per axis** — pad one axis only (the usual slab setup) instead of a
  uniform box
- **Run types** — `ENERGY_FORCE`, `GEO_OPT`, `CELL_OPT` and `MD`, with `&MOTION`
  written only when the run type needs it
- **&DFT** — MOLOPT basis sets and matching GTH potentials per `&KIND`, charge and
  multiplicity (UKS forced when needed), `&MGRID` cutoffs, `&QS EPS_DEFAULT`,
  XC functional (PBE / BLYP / BP / PADE / PBE0 / B3LYP) and DFT-D2/D3/D3(BJ)
  dispersion via `&PAIR_POTENTIAL`
- **Two SCF strategies** — `&OT` with `&OUTER_SCF` for insulators and molecules,
  or diagonalisation with Fermi-Dirac `&SMEAR`, `ADDED_MOS` and Broyden mixing for
  metals (OT cannot handle fractional occupations)
- **&KPOINTS** — Gamma-only (section omitted), a Monkhorst-Pack mesh, or a mesh
  derived from a target reciprocal-space spacing
- **&SUBSYS** — `&CELL` vectors with configurable `PERIODIC`, `&COORD` as scaled or
  cartesian coordinates, and one `&KIND` per element

The `&KIND` potential follows the functional (`GTH-PBE`, `GTH-BLYP`, `GTH-PADE`).
These are the functional aliases from `GTH_POTENTIALS`, which resolve to the right
valence for every element — no per-element `-qN` guessing.

- **Checks** — a warning strip flags the classic mistakes: a molecule sampled
  with a dense k-mesh, k-points across a slab's vacuum, too little vacuum,
  cut-offs that are too low, a relaxation with zero steps, and charged cells

- **Charge and multiplicity** can be read straight from the molecule

## Install

Plugin Manager → install from the MoleditPy plugin registry, or drop the
`cp2k_input_generator` folder into your MoleditPy plugins directory.

Requires `numpy` (already a MoleditPy dependency). `pymatgen` is optional and is
used only to expand a CIF Viewer structure that holds nothing but the asymmetric
unit; reading a `.cif` file directly never needs it.

## Use

**File → Export → CP2K Input (Quickstep)...**

Pick a structure source on the *Structure* tab, set the run up on *DFT*, *SCF* and
*Motion*, check the *Preview* tab, then **Save Input...**.

## Shared modules

`cell_model.py`, `elements.py` and `structure_panel.py` are shared byte-for-byte
with the VASP and Quantum ESPRESSO input generator plugins. Each carries its own
`SHARED_MODULE_NAME` / `SHARED_MODULE_VERSION`, independent of `PLUGIN_VERSION`:
change one, bump its version, copy it to the sibling plugins, and update the pin
in each `tests/test_structure_panel.py`.

The CIF reader, lattice construction and symmetry de-duplication are derived from
the [MoleditPy CIF Viewer](https://github.com/HiroYokoyama/moleditpy_cif_viewer)
plugin's parser.

## Tests

```bash
python -m pytest tests/ -v
```

## Licence

GNU General Public License v3.0 — see [LICENSE](LICENSE).
