import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cp2k_input_generator import writer  # noqa: E402
from cp2k_input_generator.main_dialog import Cp2kInputDialog, _to_float  # noqa: E402

from test_cell_model import _FakeMol  # noqa: E402


@pytest.fixture
def dialog(qapp):
    mol = _FakeMol(["O", "H", "H"], [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
    dlg = Cp2kInputDialog(persistent_settings=writer.default_settings(), get_molecule=lambda: mol)
    yield dlg
    dlg.deleteLater()


def test_dialog_builds_preview(dialog):
    text = dialog.preview.toPlainText()
    assert "&GLOBAL" in text and "&FORCE_EVAL" in text and "&SUBSYS" in text
    assert dialog._cell is not None
    assert dialog.save_button.isEnabled()


def test_dialog_preview_tracks_run_type(dialog):
    dialog.run_type_combo.setCurrentText("GEO_OPT (relax ions)")
    text = dialog.preview.toPlainText()
    assert "RUN_TYPE GEO_OPT" in text
    assert "&MOTION" in text


def test_dialog_preview_tracks_kpoints(dialog):
    dialog.kmode_combo.setCurrentText("Monkhorst-Pack mesh")
    dialog.kmesh_spins[0].setValue(3)
    assert "MONKHORST-PACK 3" in dialog.preview.toPlainText()


def test_dialog_updates_persistent_settings(dialog):
    dialog.cutoff_spin.setValue(650.0)
    assert dialog.persistent_settings["cutoff"] == 650.0


def test_dialog_marks_the_project_modified(qapp):
    seen = []
    mol = _FakeMol(["H"], [[0.0, 0.0, 0.0]])
    dlg = Cp2kInputDialog(
        persistent_settings={}, get_molecule=lambda: mol, mark_modified=lambda: seen.append(1)
    )
    dlg.cutoff_spin.setValue(450.0)
    assert seen
    dlg.deleteLater()


def test_dialog_settings_roundtrip(dialog):
    settings = dialog.read_settings()
    settings.update(
        {
            "run_type": "MD (molecular dynamics)",
            "cutoff": 800.0,
            "basis_set": "TZV2P-MOLOPT-GTH",
            "periodic": "XY",
        }
    )
    dialog.apply_settings(settings)
    out = dialog.read_settings()
    assert out["run_type"] == "MD (molecular dynamics)"
    assert out["cutoff"] == 800.0
    assert out["basis_set"] == "TZV2P-MOLOPT-GTH"
    assert out["periodic"] == "XY"


def test_dialog_solver_switch_toggles_fields(dialog):
    dialog.solver_combo.setCurrentText(writer.SCF_SOLVERS[1])
    assert not dialog.outer_scf_check.isEnabled()
    assert dialog.added_mos_spin.isEnabled()
    assert dialog.smear_temp_spin.isEnabled()
    assert "&SMEAR" in dialog.preview.toPlainText()

    dialog.solver_combo.setCurrentText(writer.SCF_SOLVERS[0])
    assert dialog.outer_scf_check.isEnabled()
    assert not dialog.added_mos_spin.isEnabled()
    assert "&OT" in dialog.preview.toPlainText()


def test_dialog_kind_summary_label(dialog):
    assert "O: DZVP-MOLOPT-SR-GTH / GTH-PBE" in dialog.kind_label.text()


def test_dialog_kind_summary_follows_functional(dialog):
    dialog.functional_combo.setCurrentText("BLYP")
    assert "GTH-BLYP" in dialog.kind_label.text()


def test_dialog_reports_a_missing_molecule(qapp):
    dlg = Cp2kInputDialog(persistent_settings={}, get_molecule=lambda: None)
    assert "No molecule" in dlg.preview.toPlainText()
    assert not dlg.save_button.isEnabled()
    assert dlg.kind_label.text() == "-"
    dlg.deleteLater()


@pytest.mark.parametrize(
    "text,expected", [("1e-6", 1e-6), ("0.5", 0.5), ("junk", 42.0), (None, 42.0)]
)
def test_to_float(text, expected):
    assert _to_float(text, 42.0) == expected


def test_dialog_save_writes_the_input(dialog, tmp_path, monkeypatch):
    from PyQt6.QtWidgets import QFileDialog, QMessageBox

    target = tmp_path / "run.inp"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(target), ""))
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: None)
    dialog.save_input()
    text = target.read_text(encoding="utf-8")
    assert text.startswith("&GLOBAL")
    assert "\r" not in text


def test_dialog_save_is_cancellable(dialog, tmp_path, monkeypatch):
    from PyQt6.QtWidgets import QFileDialog

    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: ("", ""))
    dialog.save_input()
    assert not list(tmp_path.iterdir())


def test_dialog_save_without_structure_warns(qapp, monkeypatch):
    from PyQt6.QtWidgets import QMessageBox

    seen = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: seen.append(a))
    dlg = Cp2kInputDialog(persistent_settings={}, get_molecule=lambda: None)
    dlg.save_input()
    assert seen
    dlg.deleteLater()


def test_dialog_copy_preview(dialog, qapp):
    dialog.copy_preview()
    assert "&GLOBAL" in qapp.clipboard().text()


def DIALOG_FACTORY(context):
    mol = _FakeMol(["O", "H", "H"], [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
    return Cp2kInputDialog(persistent_settings=writer.default_settings(), get_molecule=lambda: mol, context=context)


# -- drag and drop ----------------------------------------------------------


def test_the_dialog_accepts_a_dropped_cif(dialog, tmp_path):
    from test_structure_panel import _FakeDropEvent, _FakeMime

    path = tmp_path / "dropped.cif"
    path.write_text("data_x", encoding="utf-8")
    assert dialog.acceptDrops()
    event = _FakeDropEvent(_FakeMime([str(path)]))
    dialog.dropEvent(event)
    assert event.accepted
    assert dialog.structure_panel.cif_edit.text() == str(path)


def test_a_drag_without_a_cif_is_refused_by_the_dialog(dialog):
    from test_structure_panel import _FakeDropEvent, _FakeMime

    event = _FakeDropEvent(_FakeMime(["/tmp/notes.txt"]))
    dialog.dragEnterEvent(event)
    assert event.ignored and not event.accepted
    dialog.dropEvent(event)
    assert not event.accepted


def test_a_drag_move_follows_the_same_rule(dialog):
    from test_structure_panel import _FakeDropEvent, _FakeMime

    event = _FakeDropEvent(_FakeMime(["/tmp/x.cif"]))
    dialog.dragMoveEvent(event)
    assert event.accepted


def test_the_box_is_drawn_as_soon_as_the_dialog_opens(qapp):
    """Opening the generator should show the cell, not an empty viewer."""
    pytest.importorskip("rdkit")
    from test_structure_panel import _RecordingContext

    context = _RecordingContext()
    dlg = DIALOG_FACTORY(context)
    assert context.current_molecule is not None
    assert len(context.plotter.lines) == 12
    dlg.deleteLater()
