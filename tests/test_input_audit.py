"""Keyword and dialog regressions for the 0.9.0 input audit."""

import ast
from pathlib import Path
from unittest.mock import Mock

import pytest

from cp2k_input_generator import cell_model as cm, writer
from cp2k_input_generator.main_dialog import Cp2kInputDialog, QFileDialog, QMessageBox
from shared_fixtures import CUBIC_CIF, _FakeMol


@pytest.mark.parametrize("periodic", writer.PERIODICITY)
def test_electrostatics_follow_cell_periodicity(periodic):
    text = writer.build_input(cm.parse_cif(CUBIC_CIF), {"periodic": periodic})
    if periodic == "XYZ":
        assert "&POISSON" not in text
        return
    block = text.split("&POISSON")[1].split("&END POISSON")[0]
    assert f"PERIODIC {periodic}" in block
    assert f"PSOLVER {'MT' if periodic in ('NONE', 'XY') else 'ANALYTIC'}" in block


def test_cutoff_rounding_stays_inside_skewed_cell():
    cell = cm.parse_cif(CUBIC_CIF)
    lengths, angles = (4.0, 4.0, 10.0), (90.0, 90.0, 120.0)
    cell = cm.Cell("hex", lengths, angles, cm.cell_vectors(lengths, angles), cell.atoms)
    text = writer.build_input(cell, {"functional": "PBE0"})
    line = next(line for line in text.splitlines() if "CUTOFF_RADIUS" in line)
    radius = float(line.split()[-1])
    assert radius <= writer.hf_cutoff_radius(cell)
    assert radius == pytest.approx(writer.hf_cutoff_radius(cell), abs=1e-8)


@pytest.mark.parametrize("value", [1.23456789e-9, 0.1234567890123, 123456789.123])
def test_numeric_settings_preserve_precision(value):
    assert float(writer._format(value)) == pytest.approx(value, rel=1e-14)


def test_direct_hfx_with_mesh_is_rejected():
    with pytest.raises(ValueError, match="RI-HFX"):
        writer.build_input(cm.parse_cif(CUBIC_CIF), {
            "functional": "PBE0", "kpoint_mode": writer.KPOINT_MODES[1],
            "scf_solver": writer.SCF_SOLVERS[1],
        })


def test_zero_empty_states_with_smearing_is_rejected():
    with pytest.raises(ValueError, match="ADDED_MOS"):
        writer.build_input(cm.parse_cif(CUBIC_CIF), {
            "scf_solver": writer.SCF_SOLVERS[1], "added_mos": 0,
        })


def test_invalid_hybrid_preview_blocks_save_and_recovers(qapp, monkeypatch):
    molecule = _FakeMol(["O", "H", "H"], [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    dialog = Cp2kInputDialog(get_molecule=lambda: molecule)
    try:
        dialog.functional_combo.setCurrentText("PBE0")
        dialog.kmode_combo.setCurrentText(writer.KPOINT_MODES[1])
        assert "Invalid input" in dialog.preview.toPlainText()
        assert not dialog.save_button.isEnabled()
        file_picker = Mock()
        warning = Mock()
        monkeypatch.setattr(QFileDialog, "getSaveFileName", file_picker)
        monkeypatch.setattr(QMessageBox, "warning", warning)
        dialog.save_input()
        warning.assert_called_once()
        file_picker.assert_not_called()
        dialog.kmode_combo.setCurrentText(writer.KPOINT_MODES[0])
        assert dialog.save_button.isEnabled()
        assert "&HF" in dialog.preview.toPlainText()
    finally:
        dialog.close()


@pytest.mark.parametrize("filename", ["__init__.py", "main_dialog.py", "writer.py"])
def test_plugin_main_imports_are_module_scope(filename):
    path = Path(__file__).resolve().parents[1] / "cp2k_input_generator" / filename
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            assert not any(isinstance(child, (ast.Import, ast.ImportFrom)) for child in ast.walk(node))
