"""Submit to Cluster: the handoff to the Job Manager plugin."""

import os
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PySide6.QtWidgets import QFileDialog, QMessageBox, QWidget  # noqa: E402

from cp2k_input_generator import cluster_link, writer  # noqa: E402
from cp2k_input_generator.main_dialog import Cp2kInputDialog  # noqa: E402

from shared_fixtures import _FakeMol  # noqa: E402


class _MainWindow(QWidget):
    pass


def _job_manager(result=True):
    calls = []

    def submit_file(paths, name=""):
        calls.append((paths, name))
        if isinstance(result, Exception):
            raise result
        return result

    return SimpleNamespace(submit_file=submit_file), calls


def _plain_window(module=None):
    """A main window stand-in for the link alone (no Qt needed)."""
    records = [] if module is None else [{"name": "Job Manager", "module": module}]
    return SimpleNamespace(plugin_manager=SimpleNamespace(plugins=records))


def _main_window(module=None):
    """A real parent widget for the dialog; needs the qapp fixture."""
    mw = _MainWindow()
    records = [] if module is None else [{"name": "Job Manager", "module": module}]
    mw.plugin_manager = SimpleNamespace(plugins=records)
    return mw


@pytest.fixture
def boxes(monkeypatch):
    seen = {"information": [], "warning": []}
    for kind in seen:
        monkeypatch.setattr(QMessageBox, kind, lambda *a, _k=kind, **k: seen[_k].append(a))
    return seen


def _dialog(mw):
    mol = _FakeMol(["O", "H", "H"], [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
    dlg = Cp2kInputDialog(
        parent=mw, persistent_settings=writer.default_settings(), get_molecule=lambda: mol
    )
    dlg.test_parent = mw  # a dropped parent would delete the dialog with it
    return dlg


def _save_to(monkeypatch, path):
    asked = []

    def ask(*a, **k):
        asked.append(path)
        return (str(path), "") if path else ("", "")

    monkeypatch.setattr(QFileDialog, "getSaveFileName", ask)
    return asked


class TestClusterLink:
    def test_found_by_name_and_submit_file(self):
        module, _ = _job_manager()
        assert cluster_link.find_job_manager(_plain_window(module)) is module

    @pytest.mark.parametrize(
        "mw",
        [None, SimpleNamespace(), SimpleNamespace(plugin_manager=SimpleNamespace(plugins=["x"]))],
    )
    def test_absent(self, mw):
        assert not cluster_link.is_available(mw)

    def test_an_old_job_manager_without_submit_file_is_absent(self):
        mw = SimpleNamespace(
            plugin_manager=SimpleNamespace(plugins=[{"name": "Job Manager", "module": object()}])
        )
        assert cluster_link.find_job_manager(mw) is None

    def test_a_list_is_passed_on_and_a_single_path_as_a_string(self, tmp_path):
        module, calls = _job_manager()
        mw = _plain_window(module)
        a, b = tmp_path / "a", tmp_path / "b"
        a.write_text("a")
        b.write_text("b")
        assert cluster_link.submit_to_cluster(mw, [str(a), str(b)], name="n")
        assert cluster_link.submit_to_cluster(mw, str(a))
        assert calls == [([str(a), str(b)], "n"), (str(a), "")]

    def test_a_missing_file_is_not_submitted(self, tmp_path):
        module, calls = _job_manager()
        assert not cluster_link.submit_to_cluster(_plain_window(module), str(tmp_path / "gone"))
        assert calls == []

    def test_an_exception_counts_as_declined(self, tmp_path):
        module, _ = _job_manager(RuntimeError("boom"))
        f = tmp_path / "a"
        f.write_text("a")
        assert not cluster_link.submit_to_cluster(_plain_window(module), str(f))


class TestSubmitButton:
    def test_hidden_without_job_manager(self, qapp):
        mw = _main_window()
        dlg = _dialog(mw)
        assert not dlg.submit_button.isVisibleTo(dlg)
        dlg.deleteLater()

    def test_first_submit_saves_quietly_hands_over_and_closes(
        self, qapp, tmp_path, monkeypatch, boxes
    ):
        module, calls = _job_manager()
        mw = _main_window(module)
        dlg = _dialog(mw)
        dlg.show()
        assert dlg.submit_button.isVisibleTo(dlg)
        target = tmp_path / "water.inp"
        _save_to(monkeypatch, target)
        dlg.submit_button.click()
        assert calls == [(str(target), "water")]
        assert target.read_text(encoding="utf-8") == dlg.preview.toPlainText()
        assert boxes["information"] == [] and boxes["warning"] == []
        assert not dlg.isVisible()
        dlg.deleteLater()

    def test_an_unchanged_input_is_not_saved_again(self, qapp, tmp_path, monkeypatch, boxes):
        module, calls = _job_manager()
        dlg = _dialog(_main_window(module))
        target = tmp_path / "water.inp"
        asked = _save_to(monkeypatch, target)
        dlg.save_button.click()
        assert len(boxes["information"]) == 1
        dlg.submit_to_cluster()
        assert len(asked) == 1 and calls == [(str(target), "water")]
        dlg.deleteLater()

    def test_a_changed_setting_saves_again(self, qapp, tmp_path, monkeypatch, boxes):
        module, calls = _job_manager()
        dlg = _dialog(_main_window(module))
        target = tmp_path / "water.inp"
        asked = _save_to(monkeypatch, target)
        dlg.save_input()
        dlg.cutoff_spin.setValue(dlg.cutoff_spin.value() + 50)
        dlg.submit_to_cluster()
        assert len(asked) == 2 and calls
        assert target.read_text(encoding="utf-8") == dlg.preview.toPlainText()
        dlg.deleteLater()

    def test_cancelling_the_save_submits_nothing(self, qapp, tmp_path, monkeypatch, boxes):
        module, calls = _job_manager()
        dlg = _dialog(_main_window(module))
        dlg.show()
        _save_to(monkeypatch, "")
        dlg.submit_to_cluster()
        assert calls == [] and boxes["warning"] == [] and dlg.isVisible()
        dlg.close()
        dlg.deleteLater()

    def test_a_refusal_warns_and_keeps_the_window(self, qapp, tmp_path, monkeypatch, boxes):
        module, _ = _job_manager(False)
        dlg = _dialog(_main_window(module))
        dlg.show()
        _save_to(monkeypatch, tmp_path / "water.inp")
        dlg.submit_to_cluster()
        assert len(boxes["warning"]) == 1 and dlg.isVisible()
        dlg.close()
        dlg.deleteLater()
