"""CP2K (Quickstep) Input Generator plugin for MoleditPy."""

import logging
import os

PLUGIN_NAME = "CP2K Input Generator"
PLUGIN_VERSION = "0.5.1"
PLUGIN_AUTHOR = "HiroYokoyama"
PLUGIN_DESCRIPTION = (
    "Generate CP2K Quickstep inputs from the current molecule (in a vacuum box) "
    "or from a CIF crystal structure, with supercells, GTH/MOLOPT kinds and "
    "OT or smeared-diagonalisation SCF."
)
PLUGIN_CATEGORY = "Export"
PLUGIN_TAGS = ["DFT", "Generator"]
PLUGIN_DEPENDENCIES = ["numpy", "pyvista", "rdkit"]
PLUGIN_SUPPORTED_MOLEDITPY_VERSION = ">=4.0.0, <5.0.0"

SETTINGS_FILE = os.path.join(os.path.dirname(__file__), "settings.json")
WINDOW_ID = "cp2k_input_generator_dialog"

_context = None
_dialog_opened = False


def get_default_settings():
    from .writer import default_settings

    return default_settings()


current_settings = get_default_settings()


def _open_dialog(mw):
    """Open the dialog."""
    global _dialog_opened

    if _context is not None:
        mw = _context.get_main_window()

    from .main_dialog import Cp2kInputDialog

    if _context is not None:
        existing = _context.get_window(WINDOW_ID)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return

    def _get_molecule():
        try:
            if _context is not None:
                return _context.current_molecule
        except Exception as exc:  # pragma: no cover - host API guard
            logging.warning("%s: could not read the molecule: %s", PLUGIN_NAME, exc)
        return getattr(mw, "current_mol", None)

    def _mark_modified():
        if _context is not None:
            try:
                _context.mark_project_modified()
            except Exception:  # pragma: no cover - host API guard
                pass

    def _get_cif_viewer():
        from .structure_panel import find_cif_viewer_widget

        return find_cif_viewer_widget(mw)

    _dialog_opened = True
    dlg = Cp2kInputDialog(
        parent=mw,
        persistent_settings=current_settings,
        get_molecule=_get_molecule,
        mark_modified=_mark_modified,
        get_cif_viewer=_get_cif_viewer,
        context=_context,
    )
    if _context is not None:
        _context.register_window(WINDOW_ID, dlg)
    dlg.show()



def run(main_window):
    """Entry point for the host's automatic Plugins-menu item (manual 7.1).

    Safe alongside the entry registered in initialize(): that one goes to
    the Export menu, while a module exposing ``run`` is listed under Plugins, so the
    two do not collide.
    """
    _open_dialog(main_window)

def initialize(context):
    global _context
    _context = context

    def show_dialog():
        _open_dialog(context.get_main_window())

    context.add_export_action("CP2K Input (Quickstep)...", show_dialog)

    def save_state():
        if not _dialog_opened:
            return {}
        return {"settings": dict(current_settings)}

    def load_state(data):
        if not isinstance(data, dict):
            return
        saved = data.get("settings")
        if isinstance(saved, dict):
            current_settings.update(saved)
            dlg = context.get_window(WINDOW_ID)
            if dlg is not None:
                try:
                    dlg.apply_settings(current_settings)
                except Exception as exc:  # pragma: no cover - host API guard
                    logging.warning("%s: could not apply loaded state: %s", PLUGIN_NAME, exc)

    def handle_reset():
        global _dialog_opened
        dlg = context.get_window(WINDOW_ID)
        if dlg is not None and dlg.isVisible():
            # Leave an open dialog alone: the user may still be editing.
            return
        current_settings.clear()
        current_settings.update(get_default_settings())
        _dialog_opened = False

    context.register_save_handler(save_state)
    context.register_load_handler(load_state)
    context.register_document_reset_handler(handle_reset)
