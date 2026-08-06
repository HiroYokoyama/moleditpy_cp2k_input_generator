import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cp2k_input_generator import cell_model as cm  # noqa: E402
from cp2k_input_generator import writer  # noqa: E402


@pytest.fixture
def molecule_cell():
    return cm.cell_from_molecule(
        ["O", "H", "H"],
        [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        padding=6.0,
    )


@pytest.fixture
def bulk_cell():
    lengths, angles = (4.0, 4.0, 4.0), (90.0, 90.0, 90.0)
    lattice = cm.cell_vectors(lengths, angles)
    atom = cm.CellAtom("Cu1", "Cu", np.zeros(3), np.zeros(3))
    return cm.Cell("sc", lengths, angles, lattice, (atom,), source="cif")


@pytest.fixture
def slab_cell():
    """A slab as it arrives from the Slab Builder: a CIF with a vacuum layer."""
    lengths, angles = (4.0, 4.0, 20.0), (90.0, 90.0, 90.0)
    lattice = cm.cell_vectors(lengths, angles)
    atoms = tuple(
        cm.CellAtom(
            f"Cu{index + 1}",
            "Cu",
            np.array([0.0, 0.0, height / 20.0]),
            np.array([0.0, 0.0, height]),
        )
        for index, height in enumerate((0.0, 2.0, 4.0))
    )
    return cm.Cell("slab", lengths, angles, lattice, atoms, source="cif")


def joined(messages):
    return " | ".join(messages)


# -- CP2K -------------------------------------------------------------------


def test_no_warnings_for_a_sane_bulk_run(bulk_cell):
    assert writer.validate(bulk_cell, {"periodic": "XYZ"}) == []


def test_kpoints_with_ot_is_flagged(bulk_cell):
    messages = writer.validate(
        bulk_cell, {"kpoint_mode": "Monkhorst-Pack mesh", "kmesh": [4, 4, 4]}
    )
    assert "OT solver" in joined(messages)


def test_kpoints_with_diagonalisation_is_quiet(bulk_cell):
    messages = writer.validate(
        bulk_cell,
        {
            "kpoint_mode": "Monkhorst-Pack mesh",
            "kmesh": [4, 4, 4],
            "scf_solver": writer.SCF_SOLVERS[1],
        },
    )
    assert "OT solver" not in joined(messages)


def test_slab_sampled_through_the_vacuum(slab_cell):
    messages = writer.validate(
        slab_cell,
        {
            "kpoint_mode": "Monkhorst-Pack mesh",
            "kmesh": [4, 4, 4],
            "slab_kpoints_c1": False,
            "scf_solver": writer.SCF_SOLVERS[1],
        },
    )
    assert "vacuum direction" in joined(messages)


def test_slab_kpoint_rule_reaches_the_input(slab_cell):
    text = writer.build_input(
        slab_cell,
        {"kpoint_mode": "Monkhorst-Pack mesh", "kmesh": [4, 4, 4], "scf_solver": writer.SCF_SOLVERS[1]},
    )
    assert "MONKHORST-PACK 4 4 1" in text


def test_periodic_molecule_suggests_a_poisson_solver(molecule_cell):
    assert "&POISSON" in joined(writer.validate(molecule_cell, {"periodic": "XYZ"}))


def test_non_periodic_molecule_is_quiet(molecule_cell):
    assert "&POISSON" not in joined(writer.validate(molecule_cell, {"periodic": "NONE"}))


def test_slab_periodicity_hint(slab_cell):
    assert "PERIODIC XY" in joined(writer.validate(slab_cell, {"periodic": "XYZ"}))


def test_low_cutoffs(bulk_cell):
    messages = writer.validate(bulk_cell, {"cutoff": 200.0, "rel_cutoff": 20.0})
    assert "CUTOFF" in joined(messages)
    assert "REL_CUTOFF" in joined(messages)


def test_charged_periodic_cell(bulk_cell):
    assert "background charge" in joined(writer.validate(bulk_cell, {"charge": -1}))


def test_smearing_without_empty_states(bulk_cell):
    messages = writer.validate(
        bulk_cell, {"scf_solver": writer.SCF_SOLVERS[1], "added_mos": 0}
    )
    assert "ADDED_MOS" in joined(messages)


def test_multiplicity_note(bulk_cell):
    assert "MULTIPLICITY" in joined(writer.validate(bulk_cell, {"multiplicity": 3}))


# -- structure faults reach the CP2K warnings ---------------------------------


def test_a_partially_occupied_cif_is_flagged(bulk_cell):
    disordered = cm.Cell(
        bulk_cell.name,
        bulk_cell.lengths,
        bulk_cell.angles,
        bulk_cell.lattice,
        (cm.CellAtom("Fe1", "Fe", np.zeros(3), np.zeros(3), 0.5),),
        source="cif",
    )
    assert "partially occupied" in joined(writer.validate(disordered, {"kmesh": [8, 8, 8]}))


def test_a_left_handed_lattice_is_flagged():
    lattice = np.array([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, -4.0]])
    cell = cm.Cell("lh", (4.0, 4.0, 4.0), (90.0, 90.0, 90.0), lattice, ())
    assert "left-handed" in joined(writer.validate(cell, {"kmesh": [8, 8, 8]}))


def test_overlapping_atoms_are_flagged(bulk_cell):
    doubled = cm.Cell(
        bulk_cell.name,
        bulk_cell.lengths,
        bulk_cell.angles,
        bulk_cell.lattice,
        (
            cm.CellAtom("Cu1", "Cu", np.zeros(3), np.zeros(3)),
            cm.CellAtom("Cu2", "Cu", np.array([0.01, 0.0, 0.0]), np.array([0.04, 0.0, 0.0])),
        ),
        source="cif",
    )
    assert "closer than" in joined(writer.validate(doubled, {"kmesh": [8, 8, 8]}))
