"""CP2K (Quickstep) input writer."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

from .cell_model import Cell, kpoint_mesh_from_density, sorted_by_species

RUN_TYPES = (
    "ENERGY_FORCE (single point)",
    "GEO_OPT (relax ions)",
    "CELL_OPT (relax ions + cell)",
    "MD (molecular dynamics)",
)

FUNCTIONALS = ("PBE", "BLYP", "BP", "PADE (LDA)", "PBE0", "B3LYP")

BASIS_SETS = (
    "DZVP-MOLOPT-SR-GTH",
    "DZVP-MOLOPT-GTH",
    "TZVP-MOLOPT-GTH",
    "TZV2P-MOLOPT-GTH",
    "SZV-MOLOPT-SR-GTH",
)

SCF_SOLVERS = ("OT (insulator)", "Diagonalisation + smearing (metal)")

VDW_OPTIONS = ("None", "DFT-D3", "DFT-D3(BJ)", "DFT-D2")

PERIODICITY = ("XYZ", "XY", "XZ", "YZ", "X", "Y", "Z", "NONE")

KPOINT_MODES = ("Gamma point only", "Monkhorst-Pack mesh", "Automatic (spacing)")

# GTH_POTENTIALS lists every entry under a functional alias (e.g. "GTH-PBE"),
# so the alias is always valid and avoids guessing per-element valence counts.
_POTENTIAL_BY_FUNCTIONAL = {
    "PBE": "GTH-PBE",
    "BLYP": "GTH-BLYP",
    "BP": "GTH-BP",
    "PADE (LDA)": "GTH-PADE",
    "PBE0": "GTH-PBE",
    "B3LYP": "GTH-BLYP",
}

_XC_KEYWORD = {
    "PBE": "PBE",
    "BLYP": "BLYP",
    "BP": "BP",
    "PADE (LDA)": "PADE",
    "PBE0": "PBE0",
    "B3LYP": "B3LYP",
}


def default_settings() -> Dict:
    return {
        "project": "cp2k",
        "run_type": RUN_TYPES[0],
        "print_level": "MEDIUM",
        "functional": "PBE",
        "basis_set": BASIS_SETS[0],
        "basis_file": "BASIS_MOLOPT",
        "potential_file": "GTH_POTENTIALS",
        "charge": 0,
        "multiplicity": 1,
        "uks": False,
        "cutoff": 400.0,
        "rel_cutoff": 60.0,
        "eps_default": 1e-12,
        "eps_scf": 1e-6,
        "max_scf": 50,
        "scf_solver": SCF_SOLVERS[0],
        "outer_scf": True,
        "smear_temperature": 300.0,
        "added_mos": 20,
        "vdw": VDW_OPTIONS[0],
        "periodic": "XYZ",
        "scaled_coordinates": True,
        "kpoint_mode": KPOINT_MODES[0],
        "kmesh": [4, 4, 4],
        "kspacing": 0.03,
        "slab_kpoints_c1": True,
        "max_iter": 200,
        "max_force": 4.5e-4,
        "md_steps": 1000,
        "timestep": 0.5,
        "temperature": 300.0,
        "print_forces": True,
        "extra_dft": "",
        "supercell": [1, 1, 1],
        "padding": 6.0,
        "cubic_box": False,
    }


# --------------------------------------------------------------------------
# section tree
# --------------------------------------------------------------------------


class Section:
    """A CP2K ``&NAME ... &END NAME`` block."""

    def __init__(self, name: str, argument: str = ""):
        self.name = name.upper()
        self.argument = argument
        self.items: List[Union[Tuple[str, str], "Section"]] = []

    def keyword(self, key: str, value=None) -> "Section":
        self.items.append((key.upper(), "" if value is None else _format(value)))
        return self

    def section(self, name: str, argument: str = "") -> "Section":
        child = Section(name, argument)
        self.items.append(child)
        return child

    def add(self, child: "Section") -> "Section":
        self.items.append(child)
        return child

    def raw(self, text: str) -> "Section":
        for line in str(text or "").splitlines():
            if line.strip():
                self.items.append((line.strip(), ""))
        return self

    def render(self, indent: int = 0) -> List[str]:
        pad = "  " * indent
        header = f"{pad}&{self.name}"
        if self.argument:
            header += f" {self.argument}"
        lines = [header]
        for item in self.items:
            if isinstance(item, Section):
                lines.extend(item.render(indent + 1))
            else:
                key, value = item
                lines.append(f"{pad}  {key} {value}".rstrip())
        lines.append(f"{pad}&END {self.name}")
        return lines

    def __str__(self) -> str:
        return "\n".join(self.render())


def _format(value) -> str:
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, float):
        if value != 0 and (abs(value) < 1e-3 or abs(value) >= 1e5):
            return f"{value:.2E}"
        return f"{value:g}"
    if isinstance(value, (list, tuple)):
        return " ".join(_format(item) for item in value)
    return str(value)


def _run_type(settings: Dict) -> str:
    return str(settings.get("run_type", RUN_TYPES[0])).split()[0]


# --------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------


def build_global(settings: Dict) -> Section:
    section = Section("GLOBAL")
    section.keyword("PROJECT", settings.get("project", "cp2k"))
    section.keyword("RUN_TYPE", _run_type(settings))
    section.keyword("PRINT_LEVEL", settings.get("print_level", "MEDIUM"))
    return section


def build_dft(cell: Cell, settings: Dict) -> Section:
    dft = Section("DFT")
    dft.keyword("BASIS_SET_FILE_NAME", settings.get("basis_file", "BASIS_MOLOPT"))
    dft.keyword("POTENTIAL_FILE_NAME", settings.get("potential_file", "GTH_POTENTIALS"))
    dft.keyword("CHARGE", int(settings.get("charge", 0)))
    dft.keyword("MULTIPLICITY", int(settings.get("multiplicity", 1)))
    if settings.get("uks") or int(settings.get("multiplicity", 1)) > 1:
        dft.keyword("UKS", True)

    dft.section("QS").keyword("EPS_DEFAULT", float(settings.get("eps_default", 1e-12)))

    mgrid = dft.section("MGRID")
    mgrid.keyword("CUTOFF", float(settings.get("cutoff", 400.0)))
    mgrid.keyword("REL_CUTOFF", float(settings.get("rel_cutoff", 60.0)))
    mgrid.keyword("NGRIDS", 4)

    dft.add(build_scf(settings))
    dft.add(build_xc(settings))

    kpoints = build_kpoints(cell, settings)
    if kpoints is not None:
        dft.add(kpoints)

    dft.raw(settings.get("extra_dft", ""))
    return dft


def build_scf(settings: Dict) -> Section:
    scf = Section("SCF")
    scf.keyword("SCF_GUESS", "ATOMIC")
    scf.keyword("EPS_SCF", float(settings.get("eps_scf", 1e-6)))
    scf.keyword("MAX_SCF", int(settings.get("max_scf", 50)))

    if settings.get("scf_solver", SCF_SOLVERS[0]).startswith("OT"):
        ot = scf.section("OT")
        ot.keyword("MINIMIZER", "DIIS")
        ot.keyword("PRECONDITIONER", "FULL_SINGLE_INVERSE")
        if settings.get("outer_scf", True):
            outer = scf.section("OUTER_SCF")
            outer.keyword("EPS_SCF", float(settings.get("eps_scf", 1e-6)))
            outer.keyword("MAX_SCF", 10)
    else:
        # OT cannot be combined with fractional occupations, so a smeared metal
        # run uses traditional diagonalisation plus extra empty states.
        scf.keyword("ADDED_MOS", int(settings.get("added_mos", 20)))
        scf.section("DIAGONALIZATION").keyword("ALGORITHM", "STANDARD")
        smear = scf.section("SMEAR", "ON")
        smear.keyword("METHOD", "FERMI_DIRAC")
        smear.keyword("ELECTRONIC_TEMPERATURE", f"[K] {float(settings.get('smear_temperature', 300.0)):g}")
        mixing = scf.section("MIXING", ".TRUE.")
        mixing.keyword("METHOD", "BROYDEN_MIXING")
        mixing.keyword("ALPHA", 0.2)
        mixing.keyword("NBROYDEN", 8)
    return scf


def build_xc(settings: Dict) -> Section:
    functional = settings.get("functional", "PBE")
    xc = Section("XC")
    xc.section("XC_FUNCTIONAL", _XC_KEYWORD.get(functional, "PBE"))

    vdw = settings.get("vdw", "None")
    if vdw and vdw != "None":
        potential = xc.section("VDW_POTENTIAL")
        potential.keyword("POTENTIAL_TYPE", "PAIR_POTENTIAL")
        pair = potential.section("PAIR_POTENTIAL")
        if vdw == "DFT-D2":
            pair.keyword("TYPE", "DFTD2")
            pair.keyword("PARAMETER_FILE_NAME", "dftd3.dat")
        else:
            pair.keyword("TYPE", "DFTD3(BJ)" if "BJ" in vdw else "DFTD3")
            pair.keyword("PARAMETER_FILE_NAME", "dftd3.dat")
        pair.keyword("REFERENCE_FUNCTIONAL", _XC_KEYWORD.get(functional, "PBE"))
        pair.keyword("R_CUTOFF", 15.0)
    return xc


def build_kpoints(cell: Cell, settings: Dict) -> Optional[Section]:
    if settings.get("kpoint_mode", KPOINT_MODES[0]) == "Gamma point only":
        return None
    mesh = effective_mesh(cell, settings)
    if tuple(mesh) == (1, 1, 1):
        return None

    kpoints = Section("KPOINTS")
    kpoints.keyword("SCHEME", f"MONKHORST-PACK {mesh[0]} {mesh[1]} {mesh[2]}")
    kpoints.keyword("SYMMETRY", False)
    kpoints.keyword("FULL_GRID", False)
    return kpoints


def build_subsys(cell: Cell, settings: Dict) -> Section:
    subsys = Section("SUBSYS")

    cell_section = subsys.section("CELL")
    for name, vector in zip("ABC", cell.lattice):
        cell_section.keyword(name, " ".join(f"{float(value):.10f}" for value in vector))
    cell_section.keyword("PERIODIC", settings.get("periodic", "XYZ"))

    scaled = bool(settings.get("scaled_coordinates", True))
    coord = subsys.section("COORD")
    if scaled:
        coord.keyword("SCALED", True)
    atoms, counts = sorted_by_species(cell)
    for atom in atoms:
        values = atom.fract if scaled else atom.cart
        coord.keyword(
            atom.element, " ".join(f"{float(value):.10f}" for value in values)
        )

    basis = settings.get("basis_set", BASIS_SETS[0])
    potential = _POTENTIAL_BY_FUNCTIONAL.get(settings.get("functional", "PBE"), "GTH-PBE")
    for element, _ in counts:
        kind = subsys.section("KIND", element)
        kind.keyword("BASIS_SET", basis)
        kind.keyword("POTENTIAL", potential)
    return subsys


def build_force_eval(cell: Cell, settings: Dict) -> Section:
    force_eval = Section("FORCE_EVAL")
    force_eval.keyword("METHOD", "Quickstep")
    force_eval.add(build_dft(cell, settings))
    force_eval.add(build_subsys(cell, settings))
    if settings.get("print_forces", True):
        force_eval.section("PRINT").section("FORCES", "ON")
    return force_eval


def build_motion(settings: Dict) -> Optional[Section]:
    run_type = _run_type(settings)
    if run_type == "ENERGY_FORCE":
        return None

    motion = Section("MOTION")
    if run_type in ("GEO_OPT", "CELL_OPT"):
        optimiser = motion.section(run_type)
        optimiser.keyword("OPTIMIZER", "BFGS")
        optimiser.keyword("MAX_ITER", int(settings.get("max_iter", 200)))
        optimiser.keyword("MAX_FORCE", float(settings.get("max_force", 4.5e-4)))
    elif run_type == "MD":
        md = motion.section("MD")
        md.keyword("ENSEMBLE", "NVT")
        md.keyword("STEPS", int(settings.get("md_steps", 1000)))
        md.keyword("TIMESTEP", float(settings.get("timestep", 0.5)))
        md.keyword("TEMPERATURE", float(settings.get("temperature", 300.0)))
        thermostat = md.section("THERMOSTAT")
        thermostat.keyword("TYPE", "NOSE")
        thermostat.section("NOSE").keyword("TIMECON", 100.0)
    return motion


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------


def build_input(cell: Cell, settings: Optional[Dict] = None) -> str:
    settings = {**default_settings(), **(settings or {})}
    blocks: List[Section] = [build_global(settings), build_force_eval(cell, settings)]
    motion = build_motion(settings)
    if motion is not None:
        blocks.append(motion)
    return "\n".join(str(block) for block in blocks) + "\n"


def build_preview(cell: Cell, settings: Optional[Dict] = None) -> str:
    return build_input(cell, settings)


def suggested_filename(settings: Optional[Dict] = None) -> str:
    settings = {**default_settings(), **(settings or {})}
    project = str(settings.get("project") or "cp2k").strip() or "cp2k"
    return f"{project}.inp"


def kind_summary(cell: Cell, settings: Optional[Dict] = None) -> List[Tuple[str, str, str]]:
    """(element, basis set, potential) triples, for the UI."""
    settings = {**default_settings(), **(settings or {})}
    basis = settings.get("basis_set", BASIS_SETS[0])
    potential = _POTENTIAL_BY_FUNCTIONAL.get(settings.get("functional", "PBE"), "GTH-PBE")
    return [(element, basis, potential) for element, _ in sorted_by_species(cell)[1]]


# --------------------------------------------------------------------------
# validation
# --------------------------------------------------------------------------


def effective_mesh(cell: Cell, settings: Dict) -> Tuple[int, int, int]:
    mode = settings.get("kpoint_mode", KPOINT_MODES[0])
    if mode == "Gamma point only":
        return (1, 1, 1)
    if mode == "Automatic (spacing)":
        mesh = list(kpoint_mesh_from_density(cell, float(settings.get("kspacing", 0.03))))
    else:
        mesh = [max(1, int(value)) for value in settings.get("kmesh", [4, 4, 4])]
    if cell.source == "slab" and settings.get("slab_kpoints_c1", True):
        mesh[2] = 1
    return tuple(mesh)


def validate(cell: Cell, settings: Optional[Dict] = None) -> List[str]:
    """Warnings about setups that run but give wrong or wasteful answers."""
    settings = {**default_settings(), **(settings or {})}
    messages: List[str] = []

    mesh = effective_mesh(cell, settings)
    uses_kpoints = max(mesh) > 1
    is_ot = str(settings.get("scf_solver", SCF_SOLVERS[0])).startswith("OT")

    if uses_kpoints and is_ot:
        messages.append(
            "k-points combined with the OT solver is not supported by every CP2K build — "
            "switch to the diagonalisation solver if the run stops at the &KPOINTS section."
        )
    if cell.source == "slab" and mesh[2] > 1:
        messages.append(
            "The slab is sampled along the vacuum direction; set the third k-point to 1."
        )
    if cell.source == "molecule" and settings.get("periodic", "XYZ") == "XYZ":
        messages.append(
            "An isolated molecule is being treated as fully periodic — add a &POISSON section "
            "with PERIODIC NONE and PSOLVER MT (or WAVELET) to remove the image interaction."
        )
    if cell.source == "slab" and settings.get("periodic", "XYZ") == "XYZ":
        messages.append(
            "For a charged or strongly polar slab consider PERIODIC XY with PSOLVER MT so the "
            "vacuum is not periodic."
        )

    cutoff = float(settings.get("cutoff", 400.0))
    rel_cutoff = float(settings.get("rel_cutoff", 60.0))
    if cutoff < 300.0:
        messages.append(
            f"CUTOFF = {cutoff:g} Ry is low for MOLOPT/GTH; 400-600 Ry is typical and must be "
            "converged together with REL_CUTOFF."
        )
    if rel_cutoff < 40.0:
        messages.append(f"REL_CUTOFF = {rel_cutoff:g} Ry is low; 50-70 Ry is the usual range.")

    multiplicity = int(settings.get("multiplicity", 1))
    if multiplicity > 1 and not settings.get("uks"):
        messages.append(
            f"MULTIPLICITY {multiplicity} needs an unrestricted calculation — UKS is written "
            "automatically here, but check that it is what you intended."
        )
    if int(settings.get("charge", 0)) != 0 and settings.get("periodic", "XYZ") != "NONE":
        messages.append(
            "A charged periodic cell is neutralised by a background charge; use &POISSON with "
            "PSOLVER MT for a genuinely isolated charged system."
        )

    if not is_ot and int(settings.get("added_mos", 0)) < 1:
        messages.append(
            "Smearing needs empty states — set ADDED_MOS to roughly 20% of the occupied orbitals."
        )
    return messages
