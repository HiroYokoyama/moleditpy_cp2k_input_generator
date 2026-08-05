import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cp2k_input_generator import cell_model as cm  # noqa: E402
from cp2k_input_generator import writer  # noqa: E402


@pytest.fixture
def water_cell():
    return cm.cell_from_molecule(
        ["O", "H", "H"],
        [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        padding=5.0,
        name="water",
    )


def parse_sections(text):
    """Parse a CP2K input into a section tree; also proves every &END matches."""
    root = {"name": "", "argument": "", "entries": [], "children": []}
    stack = [root]
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("&END"):
            closing = stripped.split(None, 1)[1] if " " in stripped else ""
            assert len(stack) > 1, f"unbalanced &END {closing}"
            opened = stack.pop()
            assert opened["name"] == closing, f"&{opened['name']} closed by &END {closing}"
            continue
        if stripped.startswith("&"):
            parts = stripped[1:].split(None, 1)
            node = {
                "name": parts[0],
                "argument": parts[1] if len(parts) > 1 else "",
                "entries": [],
                "children": [],
            }
            stack[-1]["children"].append(node)
            stack.append(node)
            continue
        key, _, value = stripped.partition(" ")
        stack[-1]["entries"].append((key, value.strip()))
    assert len(stack) == 1, f"unclosed section &{stack[-1]['name']}"
    return root


def find_section(node, name, argument=None):
    for child in node["children"]:
        if child["name"] == name and (argument is None or child["argument"] == argument):
            return child
        found = find_section(child, name, argument)
        if found is not None:
            return found
    return None


def section(text, name, argument=None):
    node = find_section(parse_sections(text), name, argument)
    assert node is not None, f"section &{name} {argument or ''} missing"
    return node


def keywords(text, name, argument=None):
    """Direct keywords of a section — nested sections are not merged in."""
    return dict(section(text, name, argument)["entries"])


# -- Section primitive -----------------------------------------------------


def test_section_renders_nested_blocks():
    root = writer.Section("FORCE_EVAL")
    root.keyword("METHOD", "Quickstep")
    child = root.section("DFT")
    child.keyword("CHARGE", 0)
    assert str(root).splitlines() == [
        "&FORCE_EVAL",
        "  METHOD Quickstep",
        "  &DFT",
        "    CHARGE 0",
        "  &END DFT",
        "&END FORCE_EVAL",
    ]


def test_section_argument_is_written_on_the_header():
    assert str(writer.Section("SMEAR", "ON")).splitlines()[0] == "&SMEAR ON"


def test_section_raw_lines_are_indented():
    root = writer.Section("DFT").raw("A .TRUE.\n\nB 2")
    assert str(root).splitlines()[1:3] == ["  A .TRUE.", "  B 2"]


def test_format_helper():
    assert writer._format(True) == ".TRUE."
    assert writer._format(False) == ".FALSE."
    assert writer._format(1e-12) == "1.00E-12"
    assert writer._format(0.5) == "0.5"
    assert writer._format([1, 2, 3]) == "1 2 3"


# -- &GLOBAL ---------------------------------------------------------------


def test_global_defaults(water_cell):
    entries = keywords(writer.build_input(water_cell), "GLOBAL")
    assert entries["PROJECT"] == "cp2k"
    assert entries["RUN_TYPE"] == "ENERGY_FORCE"
    assert entries["PRINT_LEVEL"] == "MEDIUM"


@pytest.mark.parametrize(
    "run_type,expected",
    [
        ("GEO_OPT (relax ions)", "GEO_OPT"),
        ("CELL_OPT (relax ions + cell)", "CELL_OPT"),
        ("MD (molecular dynamics)", "MD"),
    ],
)
def test_global_run_type(water_cell, run_type, expected):
    entries = keywords(writer.build_input(water_cell, {"run_type": run_type}), "GLOBAL")
    assert entries["RUN_TYPE"] == expected


# -- &DFT ------------------------------------------------------------------


def test_dft_files_and_charge(water_cell):
    entries = keywords(writer.build_input(water_cell), "DFT")
    assert entries["BASIS_SET_FILE_NAME"] == "BASIS_MOLOPT"
    assert entries["POTENTIAL_FILE_NAME"] == "GTH_POTENTIALS"
    assert entries["CHARGE"] == "0"
    assert entries["MULTIPLICITY"] == "1"
    assert "UKS" not in entries


def test_dft_uks_forced_by_multiplicity(water_cell):
    entries = keywords(writer.build_input(water_cell, {"multiplicity": 3}), "DFT")
    assert entries["UKS"] == ".TRUE."


def test_dft_uks_explicit(water_cell):
    entries = keywords(writer.build_input(water_cell, {"uks": True}), "DFT")
    assert entries["UKS"] == ".TRUE."


def test_mgrid_cutoffs(water_cell):
    entries = keywords(
        writer.build_input(water_cell, {"cutoff": 600.0, "rel_cutoff": 70.0}), "MGRID"
    )
    assert entries["CUTOFF"] == "600"
    assert entries["REL_CUTOFF"] == "70"
    assert entries["NGRIDS"] == "4"


def test_qs_eps_default(water_cell):
    assert keywords(writer.build_input(water_cell), "QS")["EPS_DEFAULT"] == "1.00E-12"


def test_extra_dft_lines(water_cell):
    text = writer.build_input(water_cell, {"extra_dft": "SURFACE_DIPOLE_CORRECTION .TRUE."})
    assert keywords(text, "DFT")["SURFACE_DIPOLE_CORRECTION"] == ".TRUE."


# -- &XC -------------------------------------------------------------------


@pytest.mark.parametrize(
    "functional,expected",
    [("PBE", "PBE"), ("BLYP", "BLYP"), ("PADE (LDA)", "PADE"), ("B3LYP", "B3LYP")],
)
def test_xc_functional_section(water_cell, functional, expected):
    text = writer.build_input(water_cell, {"functional": functional})
    assert f"&XC_FUNCTIONAL {expected}" in text


def test_no_vdw_section_by_default(water_cell):
    assert "&VDW_POTENTIAL" not in writer.build_input(water_cell)


@pytest.mark.parametrize(
    "option,expected", [("DFT-D3", "DFTD3"), ("DFT-D3(BJ)", "DFTD3(BJ)"), ("DFT-D2", "DFTD2")]
)
def test_vdw_pair_potential(water_cell, option, expected):
    text = writer.build_input(water_cell, {"vdw": option})
    entries = keywords(text, "PAIR_POTENTIAL")
    assert entries["TYPE"] == expected
    assert entries["PARAMETER_FILE_NAME"] == "dftd3.dat"
    assert entries["REFERENCE_FUNCTIONAL"] == "PBE"


def test_vdw_reference_follows_the_functional(water_cell):
    text = writer.build_input(water_cell, {"vdw": "DFT-D3", "functional": "BLYP"})
    assert keywords(text, "PAIR_POTENTIAL")["REFERENCE_FUNCTIONAL"] == "BLYP"


# -- &SCF ------------------------------------------------------------------


def test_scf_defaults_to_ot(water_cell):
    text = writer.build_input(water_cell)
    entries = keywords(text, "SCF")
    assert entries["SCF_GUESS"] == "ATOMIC"
    assert entries["EPS_SCF"] == "1.00E-06"
    assert entries["MAX_SCF"] == "50"
    assert "&OT" in text
    assert "&OUTER_SCF" in text
    assert "&SMEAR" not in text


def test_scf_outer_can_be_disabled(water_cell):
    assert "&OUTER_SCF" not in writer.build_input(water_cell, {"outer_scf": False})


def test_scf_metal_solver_replaces_ot_with_smearing(water_cell):
    text = writer.build_input(
        water_cell, {"scf_solver": writer.SCF_SOLVERS[1], "added_mos": 40, "smear_temperature": 500.0}
    )
    assert "&OT" not in text
    assert "&OUTER_SCF" not in text
    assert "&DIAGONALIZATION" in text
    assert keywords(text, "SCF")["ADDED_MOS"] == "40"
    smear = keywords(text, "SMEAR")
    assert smear["METHOD"] == "FERMI_DIRAC"
    assert smear["ELECTRONIC_TEMPERATURE"] == "[K] 500"
    assert keywords(text, "MIXING")["METHOD"] == "BROYDEN_MIXING"


# -- &KPOINTS --------------------------------------------------------------


def test_no_kpoints_section_for_gamma(water_cell):
    assert "&KPOINTS" not in writer.build_input(water_cell)


def test_kpoints_monkhorst_pack(water_cell):
    text = writer.build_input(
        water_cell, {"kpoint_mode": "Monkhorst-Pack mesh", "kmesh": [2, 3, 4]}
    )
    assert keywords(text, "KPOINTS")["SCHEME"] == "MONKHORST-PACK 2 3 4"


def test_kpoints_mesh_of_one_is_omitted(water_cell):
    text = writer.build_input(
        water_cell, {"kpoint_mode": "Monkhorst-Pack mesh", "kmesh": [1, 1, 1]}
    )
    assert "&KPOINTS" not in text


def test_kpoints_from_spacing(water_cell):
    text = writer.build_input(
        water_cell, {"kpoint_mode": "Automatic (spacing)", "kspacing": 0.02}
    )
    assert keywords(text, "KPOINTS")["SCHEME"].startswith("MONKHORST-PACK")


# -- &SUBSYS ---------------------------------------------------------------


def test_cell_vectors_are_written(water_cell):
    entries = keywords(writer.build_input(water_cell), "CELL")
    for name, vector in zip("ABC", water_cell.lattice):
        assert [float(v) for v in entries[name].split()] == pytest.approx(
            list(vector), abs=1e-9
        )
    assert entries["PERIODIC"] == "XYZ"


def test_periodicity_is_configurable(water_cell):
    entries = keywords(writer.build_input(water_cell, {"periodic": "NONE"}), "CELL")
    assert entries["PERIODIC"] == "NONE"


def test_coord_scaled_by_default(water_cell):
    entries = section(writer.build_input(water_cell), "COORD")["entries"]
    assert entries[0] == ("SCALED", ".TRUE.")
    assert len(entries) == 4
    for _, value in entries[1:]:
        for number in [float(v) for v in value.split()]:
            assert 0.0 <= number <= 1.0


def test_coord_cartesian(water_cell):
    text = writer.build_input(water_cell, {"scaled_coordinates": False})
    entries = section(text, "COORD")["entries"]
    assert entries[0][0] == "O"
    assert [float(v) for v in entries[0][1].split()] == pytest.approx(
        list(water_cell.atoms[0].cart), abs=1e-9
    )


def test_coord_is_species_grouped():
    cell = cm.cell_from_molecule(
        ["H", "O", "H"], [[0, 0, 0], [1.0, 0, 0], [2.0, 0, 0]], padding=3.0
    )
    entries = section(writer.build_input(cell), "COORD")["entries"]
    assert [key for key, _ in entries[1:]] == ["H", "H", "O"]


def test_kind_sections_per_element(water_cell):
    text = writer.build_input(water_cell)
    assert "&KIND O" in text
    assert "&KIND H" in text
    assert text.count("&KIND ") == 2
    entries = keywords(text, "KIND", "O")
    assert entries["BASIS_SET"] == "DZVP-MOLOPT-SR-GTH"
    assert entries["POTENTIAL"] == "GTH-PBE"


@pytest.mark.parametrize(
    "functional,potential",
    [("PBE", "GTH-PBE"), ("BLYP", "GTH-BLYP"), ("PADE (LDA)", "GTH-PADE"), ("B3LYP", "GTH-BLYP")],
)
def test_kind_potential_matches_functional(water_cell, functional, potential):
    entries = keywords(writer.build_input(water_cell, {"functional": functional}), "KIND", "O")
    assert entries["POTENTIAL"] == potential


def test_kind_summary_helper(water_cell):
    assert writer.kind_summary(water_cell) == [
        ("O", "DZVP-MOLOPT-SR-GTH", "GTH-PBE"),
        ("H", "DZVP-MOLOPT-SR-GTH", "GTH-PBE"),
    ]


# -- &MOTION ---------------------------------------------------------------


def test_no_motion_section_for_single_point(water_cell):
    assert "&MOTION" not in writer.build_input(water_cell)


def test_geo_opt_motion(water_cell):
    text = writer.build_input(
        water_cell, {"run_type": "GEO_OPT (relax ions)", "max_iter": 50, "max_force": 1e-3}
    )
    entries = keywords(text, "GEO_OPT")
    assert entries["OPTIMIZER"] == "BFGS"
    assert entries["MAX_ITER"] == "50"
    assert entries["MAX_FORCE"] == "0.001"


def test_cell_opt_motion(water_cell):
    text = writer.build_input(water_cell, {"run_type": "CELL_OPT (relax ions + cell)"})
    assert "&CELL_OPT" in text
    assert "&GEO_OPT" not in text


def test_md_motion(water_cell):
    text = writer.build_input(
        water_cell,
        {"run_type": "MD (molecular dynamics)", "md_steps": 500, "timestep": 1.0, "temperature": 450.0},
    )
    entries = keywords(text, "MD")
    assert entries["ENSEMBLE"] == "NVT"
    assert entries["STEPS"] == "500"
    assert entries["TIMESTEP"] == "1"
    assert entries["TEMPERATURE"] == "450"
    assert keywords(text, "THERMOSTAT")["TYPE"] == "NOSE"


def test_print_forces_toggle(water_cell):
    assert "&FORCES ON" in writer.build_input(water_cell)
    assert "&FORCES ON" not in writer.build_input(water_cell, {"print_forces": False})


# -- assembly --------------------------------------------------------------


def test_every_section_is_closed(water_cell):
    """parse_sections asserts that every &SECTION is closed by a matching &END."""
    text = writer.build_input(water_cell, {"run_type": "MD (molecular dynamics)", "vdw": "DFT-D3"})
    root = parse_sections(text)
    assert [child["name"] for child in root["children"]] == ["GLOBAL", "FORCE_EVAL", "MOTION"]


def test_block_order(water_cell):
    text = writer.build_input(water_cell, {"run_type": "GEO_OPT (relax ions)"})
    order = [text.index(token) for token in ("&GLOBAL", "&FORCE_EVAL", "&DFT", "&SUBSYS", "&MOTION")]
    assert order == sorted(order)


def test_input_ends_with_newline(water_cell):
    assert writer.build_input(water_cell).endswith("\n")


def test_suggested_filename():
    assert writer.suggested_filename({"project": "si_bulk"}) == "si_bulk.inp"
    assert writer.suggested_filename({"project": ""}) == "cp2k.inp"


def test_default_settings_are_independent_copies():
    first = writer.default_settings()
    first["cutoff"] = 1.0
    assert writer.default_settings()["cutoff"] == 400.0
