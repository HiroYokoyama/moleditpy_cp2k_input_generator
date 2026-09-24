"""Regressions fixed in 0.8.0; each test fails on 0.7.3."""

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cp2k_input_generator import cell_model as cm  # noqa: E402
from cp2k_input_generator import writer  # noqa: E402
from cp2k_input_generator.main_dialog import Cp2kInputDialog  # noqa: E402
from cp2k_input_generator.structure_panel import SOURCE_CIF  # noqa: E402

from shared_fixtures import CUBIC_CIF, _FakeMol  # noqa: E402


def _hexagonal(a=4.0, c=10.0):
    lengths, angles = (a, a, c), (90.0, 90.0, 120.0)
    atom = cm.CellAtom("C1", "C", [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    return cm.Cell("hex", lengths, angles, cm.cell_vectors(lengths, angles), (atom,))


def test_hf_radius_fits_the_inscribed_sphere():
    """Half the shortest face-to-face width: a*sin(60)/2, not a/2."""
    assert writer.hf_cutoff_radius(_hexagonal()) == pytest.approx(4.0 * 3 ** 0.5 / 4, abs=1e-6)


def test_hf_radius_of_a_cube_is_unchanged():
    assert writer.hf_cutoff_radius(cm.parse_cif(CUBIC_CIF)) == pytest.approx(2.0)


@pytest.mark.parametrize("periodic", ["NONE", "XY"])
def test_poisson_follows_a_non_periodic_cell(periodic):
    """&CELL PERIODIC alone leaves the electrostatics fully periodic."""
    text = writer.build_input(cm.parse_cif(CUBIC_CIF), {"periodic": periodic})
    assert "&POISSON" in text
    assert f"PERIODIC {periodic}" in text.split("&POISSON")[1].split("&END POISSON")[0]
    assert "PSOLVER MT" in text


def test_fully_periodic_cell_writes_no_poisson():
    assert "&POISSON" not in writer.build_input(cm.parse_cif(CUBIC_CIF), {"periodic": "XYZ"})


def test_unsupported_partial_periodicity_is_flagged():
    messages = writer.validate(cm.parse_cif(CUBIC_CIF), {"periodic": "X"})
    assert any("&POISSON" in message for message in messages)


class _RadicalAtom:
    def __init__(self, symbol, radicals):
        self.symbol, self.radicals = symbol, radicals

    def GetSymbol(self):
        return self.symbol

    def HasProp(self, name):
        return False

    def GetFormalCharge(self):
        return 0

    def GetNumRadicalElectrons(self):
        return self.radicals


class _Methyl(_FakeMol):
    def __init__(self):
        super().__init__(["C", "H", "H", "H"], [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0]])
        self._atoms = [_RadicalAtom("C", 1)] + [_RadicalAtom("H", 0)] * 3


@pytest.fixture
def dialog(qapp, tmp_path):
    cif = tmp_path / "bulk.cif"
    cif.write_text(CUBIC_CIF, encoding="utf-8")
    mol = _Methyl()
    dlg = Cp2kInputDialog(persistent_settings=writer.default_settings(), get_molecule=lambda: mol)
    dlg.auto_charge_check.setChecked(True)
    yield dlg, str(cif)
    dlg.deleteLater()


def test_auto_multiplicity_reads_the_boxed_molecule(dialog):
    dlg, _ = dialog
    dlg.update_preview()
    assert dlg.multiplicity_spin.value() == 2


def test_auto_multiplicity_ignores_the_molecule_for_a_cif(dialog):
    """A doublet radical in the editor must not make a CIF crystal a doublet."""
    dlg, path = dialog
    dlg.multiplicity_spin.setValue(1)
    dlg.structure_panel.source_combo.setCurrentText(SOURCE_CIF)
    dlg.structure_panel.cif_edit.setText(path)
    dlg.update_preview()
    assert dlg.multiplicity_spin.value() == 1
    assert "UKS" not in dlg.preview.toPlainText()
