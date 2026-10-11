"""CP2K input generator dialog."""

from __future__ import annotations

import os

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QSpinBox,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from . import PLUGIN_NAME, PLUGIN_VERSION
from . import writer
from .cell_model import molecule_charge_and_multiplicity
from .structure_panel import SOURCE_MOLECULE, StructurePanel, dropped_cif_path


class Cp2kInputDialog(QDialog):
    def __init__(
        self,
        parent=None,
        persistent_settings=None,
        get_molecule=None,
        mark_modified=None,
        get_cif_viewer=None,
        context=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(f"{PLUGIN_NAME} v{PLUGIN_VERSION}")
        self.resize(940, 720)

        self.persistent_settings = persistent_settings if persistent_settings is not None else {}
        self.mark_modified = mark_modified
        self.context = context
        self.setAcceptDrops(True)
        self._updating = False
        self._cell = None
        self._get_molecule = get_molecule

        layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, 1)

        self.structure_panel = StructurePanel(
            get_molecule=get_molecule, get_cif_viewer=get_cif_viewer, context=context
        )
        structure_tab = QWidget()
        structure_layout = QVBoxLayout(structure_tab)
        structure_layout.addWidget(self.structure_panel)
        structure_layout.addWidget(self._build_subsys_box())
        structure_layout.addStretch(1)
        self.tabs.addTab(structure_tab, "Structure")

        self.tabs.addTab(self._build_dft_tab(), "DFT")
        self.tabs.addTab(self._build_scf_tab(), "SCF")
        self.tabs.addTab(self._build_motion_tab(), "Motion")

        self.preview = QPlainTextEdit()
        self.preview.setReadOnly(True)
        self.preview.setFont(QFont("Courier New", 9))
        self.tabs.addTab(self.preview, "Preview")

        self.warning_label = QLabel()
        self.warning_label.setWordWrap(True)
        self.warning_label.setTextFormat(Qt.TextFormat.RichText)
        self.warning_label.setStyleSheet("QLabel { color: #b36b00; }")
        self.warning_label.setVisible(False)
        layout.addWidget(self.warning_label)

        buttons = QDialogButtonBox()
        self.save_button = buttons.addButton("Save Input...", QDialogButtonBox.ButtonRole.AcceptRole)
        self.copy_button = buttons.addButton("Copy", QDialogButtonBox.ButtonRole.ActionRole)
        buttons.addButton(QDialogButtonBox.StandardButton.Close)
        self.save_button.clicked.connect(self.save_input)
        self.copy_button.clicked.connect(self.copy_preview)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

        self.apply_settings(self.persistent_settings)
        self.structure_panel.changed.connect(self.update_preview)
        self.update_preview()

    # -- tabs -------------------------------------------------------------

    def _build_subsys_box(self) -> QGroupBox:
        box = QGroupBox("&SUBSYS")
        form = QFormLayout(box)
        self.periodic_combo = QComboBox()
        self.periodic_combo.addItems(writer.PERIODICITY)
        form.addRow("PERIODIC:", self.periodic_combo)
        self.scaled_check = QCheckBox("Write &COORD as SCALED (fractional) coordinates")
        form.addRow("", self.scaled_check)
        self.kind_label = QLabel("-")
        self.kind_label.setWordWrap(True)
        form.addRow("&KIND:", self.kind_label)

        self.periodic_combo.currentTextChanged.connect(self.update_preview)
        self.scaled_check.toggled.connect(self.update_preview)
        return box

    def _build_dft_tab(self) -> QWidget:
        tab = QWidget()
        outer = QVBoxLayout(tab)

        box = QGroupBox("&GLOBAL")
        form = QFormLayout(box)
        self.project_edit = QLineEdit()
        form.addRow("PROJECT:", self.project_edit)
        self.run_type_combo = QComboBox()
        self.run_type_combo.addItems(writer.RUN_TYPES)
        form.addRow("RUN_TYPE:", self.run_type_combo)
        self.print_level_combo = QComboBox()
        self.print_level_combo.addItems(["SILENT", "LOW", "MEDIUM", "HIGH", "DEBUG"])
        form.addRow("PRINT_LEVEL:", self.print_level_combo)
        outer.addWidget(box)

        dft_box = QGroupBox("&DFT")
        dft_form = QFormLayout(dft_box)
        self.functional_combo = QComboBox()
        self.functional_combo.addItems(writer.FUNCTIONALS)
        dft_form.addRow("XC_FUNCTIONAL:", self.functional_combo)
        self.basis_combo = QComboBox()
        self.basis_combo.addItems(writer.BASIS_SETS)
        dft_form.addRow("BASIS_SET:", self.basis_combo)
        self.basis_file_edit = QLineEdit()
        dft_form.addRow("BASIS_SET_FILE_NAME:", self.basis_file_edit)
        self.potential_file_edit = QLineEdit()
        dft_form.addRow("POTENTIAL_FILE_NAME:", self.potential_file_edit)
        self.charge_spin = QSpinBox()
        self.charge_spin.setRange(-20, 20)
        dft_form.addRow("CHARGE:", self.charge_spin)
        self.multiplicity_spin = QSpinBox()
        self.multiplicity_spin.setRange(1, 20)
        dft_form.addRow("MULTIPLICITY:", self.multiplicity_spin)
        self.uks_check = QCheckBox("UKS (unrestricted)")
        dft_form.addRow("", self.uks_check)
        self.auto_charge_check = QCheckBox("Read charge and multiplicity from the molecule")
        dft_form.addRow("", self.auto_charge_check)
        self.cutoff_spin = QDoubleSpinBox()
        self.cutoff_spin.setRange(50.0, 5000.0)
        self.cutoff_spin.setSingleStep(50.0)
        self.cutoff_spin.setSuffix(" Ry")
        dft_form.addRow("MGRID CUTOFF:", self.cutoff_spin)
        self.rel_cutoff_spin = QDoubleSpinBox()
        self.rel_cutoff_spin.setRange(10.0, 500.0)
        self.rel_cutoff_spin.setSingleStep(5.0)
        self.rel_cutoff_spin.setSuffix(" Ry")
        dft_form.addRow("REL_CUTOFF:", self.rel_cutoff_spin)
        self.vdw_combo = QComboBox()
        self.vdw_combo.addItems(writer.VDW_OPTIONS)
        dft_form.addRow("Dispersion:", self.vdw_combo)
        outer.addWidget(dft_box)

        kpoint_box = QGroupBox("&KPOINTS")
        kpoint_form = QFormLayout(kpoint_box)
        self.kmode_combo = QComboBox()
        self.kmode_combo.addItems(writer.KPOINT_MODES)
        kpoint_form.addRow("Mode:", self.kmode_combo)
        mesh_widget = QWidget()
        grid = QGridLayout(mesh_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        self.kmesh_spins = []
        for column, axis in enumerate("123"):
            spin = QSpinBox()
            spin.setRange(1, 200)
            grid.addWidget(QLabel(f"n{axis}:"), 0, column * 2)
            grid.addWidget(spin, 0, column * 2 + 1)
            self.kmesh_spins.append(spin)
        kpoint_form.addRow("Mesh:", mesh_widget)
        self.kspacing_spin = QDoubleSpinBox()
        self.kspacing_spin.setRange(0.001, 1.0)
        self.kspacing_spin.setDecimals(4)
        self.kspacing_spin.setSingleStep(0.005)
        self.kspacing_spin.setSuffix(" 1/A (no 2\u03c0)")
        self.kspacing_spin.setToolTip(
            "Largest distance between k-points along each reciprocal lattice vector, "
            "without the 2\u03c0 factor: n = ceil(|b| / spacing).\n"
            "VASP's KSPACING and pymatgen include 2\u03c0, so their 0.2 is about 0.032 here."
        )
        kpoint_form.addRow("Automatic spacing:", self.kspacing_spin)
        outer.addWidget(kpoint_box)

        extra_box = QGroupBox("Additional &DFT lines")
        extra_layout = QVBoxLayout(extra_box)
        self.extra_edit = QPlainTextEdit()
        self.extra_edit.setPlaceholderText("SURFACE_DIPOLE_CORRECTION .TRUE.")
        self.extra_edit.setMaximumHeight(90)
        extra_layout.addWidget(self.extra_edit)
        outer.addWidget(extra_box)
        outer.addStretch(1)

        for widget in (
            self.run_type_combo,
            self.print_level_combo,
            self.functional_combo,
            self.basis_combo,
            self.vdw_combo,
            self.kmode_combo,
        ):
            widget.currentTextChanged.connect(self.update_preview)
        for widget in (self.project_edit, self.basis_file_edit, self.potential_file_edit):
            widget.textChanged.connect(self.update_preview)
        for widget in (self.charge_spin, self.multiplicity_spin, *self.kmesh_spins):
            widget.valueChanged.connect(self.update_preview)
        for widget in (self.cutoff_spin, self.rel_cutoff_spin, self.kspacing_spin):
            widget.valueChanged.connect(self.update_preview)
        self.uks_check.toggled.connect(self.update_preview)
        self.auto_charge_check.toggled.connect(self.update_preview)
        self.extra_edit.textChanged.connect(self.update_preview)
        return tab

    def _build_scf_tab(self) -> QWidget:
        tab = QWidget()
        outer = QVBoxLayout(tab)
        box = QGroupBox("&SCF")
        form = QFormLayout(box)
        self.solver_combo = QComboBox()
        self.solver_combo.addItems(writer.SCF_SOLVERS)
        form.addRow("Solver:", self.solver_combo)
        self.eps_scf_edit = QLineEdit()
        form.addRow("EPS_SCF:", self.eps_scf_edit)
        self.max_scf_spin = QSpinBox()
        self.max_scf_spin.setRange(1, 10000)
        form.addRow("MAX_SCF:", self.max_scf_spin)
        self.outer_scf_check = QCheckBox("&OUTER_SCF (OT only)")
        form.addRow("", self.outer_scf_check)
        self.eps_default_edit = QLineEdit()
        form.addRow("QS EPS_DEFAULT:", self.eps_default_edit)
        self.added_mos_spin = QSpinBox()
        self.added_mos_spin.setRange(0, 5000)
        form.addRow("ADDED_MOS:", self.added_mos_spin)
        self.smear_temp_spin = QDoubleSpinBox()
        self.smear_temp_spin.setRange(1.0, 100000.0)
        self.smear_temp_spin.setSuffix(" K")
        form.addRow("Smearing temperature:", self.smear_temp_spin)
        outer.addWidget(box)

        hint = QLabel(
            "OT is the fast default for insulators and molecules. Metals need "
            "fractional occupations, which OT cannot do — the second solver "
            "switches to diagonalisation with Fermi-Dirac smearing, extra empty "
            "states and Broyden mixing."
        )
        hint.setWordWrap(True)
        outer.addWidget(hint)
        outer.addStretch(1)

        self.solver_combo.currentTextChanged.connect(self.update_preview)
        for widget in (self.eps_scf_edit, self.eps_default_edit):
            widget.textChanged.connect(self.update_preview)
        for widget in (self.max_scf_spin, self.added_mos_spin):
            widget.valueChanged.connect(self.update_preview)
        self.smear_temp_spin.valueChanged.connect(self.update_preview)
        self.outer_scf_check.toggled.connect(self.update_preview)
        return tab

    def _build_motion_tab(self) -> QWidget:
        tab = QWidget()
        outer = QVBoxLayout(tab)

        opt_box = QGroupBox("&GEO_OPT / &CELL_OPT")
        opt_form = QFormLayout(opt_box)
        self.max_iter_spin = QSpinBox()
        self.max_iter_spin.setRange(1, 100000)
        opt_form.addRow("MAX_ITER:", self.max_iter_spin)
        self.max_force_edit = QLineEdit()
        opt_form.addRow("MAX_FORCE:", self.max_force_edit)
        outer.addWidget(opt_box)

        md_box = QGroupBox("&MD")
        md_form = QFormLayout(md_box)
        self.md_steps_spin = QSpinBox()
        self.md_steps_spin.setRange(1, 10000000)
        md_form.addRow("STEPS:", self.md_steps_spin)
        self.timestep_spin = QDoubleSpinBox()
        self.timestep_spin.setRange(0.01, 100.0)
        self.timestep_spin.setSingleStep(0.1)
        self.timestep_spin.setSuffix(" fs")
        md_form.addRow("TIMESTEP:", self.timestep_spin)
        self.temperature_spin = QDoubleSpinBox()
        self.temperature_spin.setRange(0.0, 100000.0)
        self.temperature_spin.setSuffix(" K")
        md_form.addRow("TEMPERATURE:", self.temperature_spin)
        outer.addWidget(md_box)

        self.print_forces_check = QCheckBox("Print forces (&FORCES ON)")
        outer.addWidget(self.print_forces_check)
        outer.addStretch(1)

        self.max_iter_spin.valueChanged.connect(self.update_preview)
        self.max_force_edit.textChanged.connect(self.update_preview)
        self.md_steps_spin.valueChanged.connect(self.update_preview)
        for widget in (self.timestep_spin, self.temperature_spin):
            widget.valueChanged.connect(self.update_preview)
        self.print_forces_check.toggled.connect(self.update_preview)
        return tab

    # -- drag and drop ----------------------------------------------------

    def dragEnterEvent(self, event) -> None:  # noqa: N802 - Qt naming
        """Accept a CIF dropped anywhere on the dialog, not only on the panel."""
        if dropped_cif_path(event.mimeData()):
            event.acceptProposedAction()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802 - Qt naming
        self.dragEnterEvent(event)

    def dropEvent(self, event) -> None:  # noqa: N802 - Qt naming
        path = dropped_cif_path(event.mimeData())
        if not path:
            event.ignore()
            return
        self.structure_panel.load_cif_path(path)
        event.acceptProposedAction()

    # -- settings ---------------------------------------------------------

    def apply_settings(self, settings) -> None:
        settings = {**writer.default_settings(), **(settings or {})}
        self._updating = True
        try:
            self.project_edit.setText(str(settings.get("project", "cp2k")))
            self.run_type_combo.setCurrentText(settings.get("run_type", writer.RUN_TYPES[0]))
            self.print_level_combo.setCurrentText(settings.get("print_level", "MEDIUM"))
            self.functional_combo.setCurrentText(settings.get("functional", "PBE"))
            self.basis_combo.setCurrentText(settings.get("basis_set", writer.BASIS_SETS[0]))
            self.basis_file_edit.setText(str(settings.get("basis_file", "BASIS_MOLOPT")))
            self.potential_file_edit.setText(str(settings.get("potential_file", "GTH_POTENTIALS")))
            self.charge_spin.setValue(int(settings.get("charge", 0)))
            self.multiplicity_spin.setValue(int(settings.get("multiplicity", 1)))
            self.uks_check.setChecked(bool(settings.get("uks")))
            self.auto_charge_check.setChecked(bool(settings.get("auto_charge", False)))
            self.cutoff_spin.setValue(float(settings.get("cutoff", 400.0)))
            self.rel_cutoff_spin.setValue(float(settings.get("rel_cutoff", 60.0)))
            self.vdw_combo.setCurrentText(settings.get("vdw", writer.VDW_OPTIONS[0]))
            self.kmode_combo.setCurrentText(settings.get("kpoint_mode", writer.KPOINT_MODES[0]))
            for spin, value in zip(self.kmesh_spins, settings.get("kmesh", [4, 4, 4])):
                spin.setValue(max(1, int(value)))
            self.kspacing_spin.setValue(float(settings.get("kspacing", 0.03)))
            self.extra_edit.setPlainText(str(settings.get("extra_dft", "") or ""))

            self.solver_combo.setCurrentText(settings.get("scf_solver", writer.SCF_SOLVERS[0]))
            self.eps_scf_edit.setText(str(settings.get("eps_scf", 1e-6)))
            self.max_scf_spin.setValue(int(settings.get("max_scf", 50)))
            self.outer_scf_check.setChecked(bool(settings.get("outer_scf", True)))
            self.eps_default_edit.setText(str(settings.get("eps_default", 1e-12)))
            self.added_mos_spin.setValue(int(settings.get("added_mos", 20)))
            self.smear_temp_spin.setValue(float(settings.get("smear_temperature", 300.0)))

            self.max_iter_spin.setValue(int(settings.get("max_iter", 200)))
            self.max_force_edit.setText(str(settings.get("max_force", 4.5e-4)))
            self.md_steps_spin.setValue(int(settings.get("md_steps", 1000)))
            self.timestep_spin.setValue(float(settings.get("timestep", 0.5)))
            self.temperature_spin.setValue(float(settings.get("temperature", 300.0)))
            self.print_forces_check.setChecked(bool(settings.get("print_forces", True)))

            self.periodic_combo.setCurrentText(settings.get("periodic", "XYZ"))
            self.scaled_check.setChecked(bool(settings.get("scaled_coordinates", True)))
            self.structure_panel.apply_settings(settings)
        finally:
            self._updating = False
        self.update_preview()

    def read_settings(self) -> dict:
        settings = {
            "project": self.project_edit.text(),
            "run_type": self.run_type_combo.currentText(),
            "print_level": self.print_level_combo.currentText(),
            "functional": self.functional_combo.currentText(),
            "basis_set": self.basis_combo.currentText(),
            "basis_file": self.basis_file_edit.text(),
            "potential_file": self.potential_file_edit.text(),
            "charge": self.charge_spin.value(),
            "multiplicity": self.multiplicity_spin.value(),
            "uks": self.uks_check.isChecked(),
            "auto_charge": self.auto_charge_check.isChecked(),
            "cutoff": self.cutoff_spin.value(),
            "rel_cutoff": self.rel_cutoff_spin.value(),
            "vdw": self.vdw_combo.currentText(),
            "kpoint_mode": self.kmode_combo.currentText(),
            "kmesh": [spin.value() for spin in self.kmesh_spins],
            "kspacing": self.kspacing_spin.value(),
            "extra_dft": self.extra_edit.toPlainText(),
            "scf_solver": self.solver_combo.currentText(),
            "eps_scf": _to_float(self.eps_scf_edit.text(), 1e-6),
            "max_scf": self.max_scf_spin.value(),
            "outer_scf": self.outer_scf_check.isChecked(),
            "eps_default": _to_float(self.eps_default_edit.text(), 1e-12),
            "added_mos": self.added_mos_spin.value(),
            "smear_temperature": self.smear_temp_spin.value(),
            "max_iter": self.max_iter_spin.value(),
            "max_force": _to_float(self.max_force_edit.text(), 4.5e-4),
            "md_steps": self.md_steps_spin.value(),
            "timestep": self.timestep_spin.value(),
            "temperature": self.temperature_spin.value(),
            "print_forces": self.print_forces_check.isChecked(),
            "periodic": self.periodic_combo.currentText(),
            "scaled_coordinates": self.scaled_check.isChecked(),
        }
        settings.update(self.structure_panel.read_settings())
        return settings

    # -- preview / output -------------------------------------------------

    def update_preview(self, *_args) -> None:
        if self._updating:
            return
        self._apply_auto_charge()
        settings = self.read_settings()
        self.persistent_settings.update(settings)
        if self.mark_modified is not None:
            try:
                self.mark_modified()
            except Exception:  # pragma: no cover - host API guard
                pass

        is_ot = self.solver_combo.currentText().startswith("OT")
        self.outer_scf_check.setEnabled(is_ot)
        self.added_mos_spin.setEnabled(not is_ot)
        self.smear_temp_spin.setEnabled(not is_ot)

        try:
            self._cell = self.structure_panel.build_cell()
        except (ValueError, OSError) as exc:
            self._cell = None
            self.structure_panel.refresh_summary(error=str(exc))
            self.kind_label.setText("-")
            self.preview.setPlainText(f"# {exc}")
            self._show_warnings([])
            self.save_button.setEnabled(False)
            return

        self.structure_panel.refresh_summary(self._cell)
        self.kind_label.setText(
            "<br>".join(
                f"{element}: {basis} / {potential}"
                for element, basis, potential in writer.kind_summary(self._cell, settings)
            )
        )
        try:
            content = writer.build_input(self._cell, settings)
        except ValueError as exc:
            self.preview.setPlainText(f"# Invalid input: {exc}")
            self._show_warnings(writer.validate(self._cell, settings))
            self.save_button.setEnabled(False)
            return
        self.preview.setPlainText(content)
        self._show_warnings(writer.validate(self._cell, settings))
        self.save_button.setEnabled(True)

    def _apply_auto_charge(self) -> None:
        if not self.auto_charge_check.isChecked() or self._get_molecule is None:
            return
        # The molecule open in MoleditPy only describes a cell boxed from it.
        # For a CIF its charge and spin belong to something else
        # entirely, so a crystal is written neutral rather than keeping
        # whatever the molecule last put in the box.
        if self.structure_panel.source_combo.currentText() != SOURCE_MOLECULE:
            charge, multiplicity = 0, 1
        else:
            try:
                charge, multiplicity = molecule_charge_and_multiplicity(self._get_molecule())
            except (ValueError, AttributeError, TypeError):
                return
        blocked = self._updating
        self._updating = True
        try:
            self.charge_spin.setValue(int(charge))
            self.multiplicity_spin.setValue(int(multiplicity))
        finally:
            self._updating = blocked

    def _show_warnings(self, messages) -> None:
        if not messages:
            self.warning_label.setVisible(False)
            self.warning_label.clear()
            return
        items = "".join(f"<li>{message}</li>" for message in messages)
        self.warning_label.setText(f"<b>Check:</b><ul>{items}</ul>")
        self.warning_label.setVisible(True)

    def copy_preview(self) -> None:
        clipboard = QApplication.clipboard()
        if clipboard is not None:
            clipboard.setText(self.preview.toPlainText())

    def save_input(self) -> None:
        if self._cell is None:
            QMessageBox.warning(self, "CP2K Input Generator", "There is no valid structure to write.")
            return
        settings = self.read_settings()
        try:
            content = writer.build_input(self._cell, settings)
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid CP2K input", str(exc))
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save CP2K input", writer.suggested_filename(settings), "CP2K input (*.inp);;All files (*)"
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(content)
        except OSError as exc:
            QMessageBox.critical(self, "CP2K Input Generator", f"Could not write the file:\n{exc}")
            return
        QMessageBox.information(self, "CP2K Input Generator", f"Wrote\n{os.path.abspath(path)}")


def _to_float(text, default):
    try:
        return float(str(text))
    except (TypeError, ValueError):
        return default
