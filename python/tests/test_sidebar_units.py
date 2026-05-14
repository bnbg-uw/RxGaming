from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import unittest
import importlib.util

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

try:
    from PySide6.QtCore import Qt  # type: ignore
    from PySide6.QtWidgets import QApplication  # type: ignore
except ModuleNotFoundError:
    Qt = None
    QApplication = None

if QApplication is not None:
    from gaming_ui.sidebar import StructureInfo, UnitListModel  # noqa: E402
    units_spec = importlib.util.spec_from_file_location("gaming_ui.units", ROOT / "python" / "gaming_ui" / "units.py")
    if units_spec is None or units_spec.loader is None:
        raise RuntimeError("Could not load gaming_ui.units")
    units = importlib.util.module_from_spec(units_spec)
    units_spec.loader.exec_module(units)


@dataclass
class FakeStructure:
    tph: float
    ba: float
    mcs: float
    cc: float


class FakeUnit:
    def __init__(self) -> None:
        self.name = "Unit 1"
        self.areaHa = 4046.8564224 / 10000.0
        self.result = SimpleNamespace(name="success")
        self.currentStructure = FakeStructure(247.105381, 30.0, 4.0, 0.5)
        self.targetStructure = FakeStructure(247.105381, 30.0, 4.0, 0.5)
        self.treatedStructure = None


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestSidebarUnits(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])

    def test_target_tph_edit_is_converted_back_to_metric(self) -> None:
        widget = StructureInfo(units.UnitSystem.IMPERIAL)
        unit = FakeUnit()
        widget.update_for_unit(unit)

        widget.target_tph.setText("100.00")
        widget._apply_target_change(0)

        self.assertAlmostEqual(100.0 / units.TPA_PER_TPH, unit.targetStructure.tph)

    def test_target_ba_edit_is_converted_back_to_metric(self) -> None:
        widget = StructureInfo(units.UnitSystem.IMPERIAL)
        unit = FakeUnit()
        widget.update_for_unit(unit)

        widget.target_ba.setText("50.00")
        widget._apply_target_change(1)

        self.assertAlmostEqual(50.0 / units.FT2_PER_ACRE_PER_M2_PER_HA, unit.targetStructure.ba)

    def test_target_cc_edit_is_converted_back_to_canonical_proportion(self) -> None:
        widget = StructureInfo(units.UnitSystem.IMPERIAL)
        unit = FakeUnit()
        widget.update_for_unit(unit)

        widget.target_cc.setText("65.00")
        widget._apply_target_change(3)

        self.assertAlmostEqual(0.65, unit.targetStructure.cc)

    def test_current_cc_displays_as_percent(self) -> None:
        widget = StructureInfo(units.UnitSystem.IMPERIAL)
        unit = FakeUnit()
        widget.update_for_unit(unit)

        self.assertEqual("50.00", widget.current_cc.text())

    def test_unit_tooltip_includes_stand_area(self) -> None:
        unit = FakeUnit()
        model = UnitListModel([unit], units.UnitSystem.IMPERIAL)

        tooltip = model.data(model.index(0, 0), role=Qt.ItemDataRole.ToolTipRole)

        self.assertIn("Area: 1.00 ac", tooltip)
        self.assertIn("BA: 130.68", tooltip)


if __name__ == "__main__":
    unittest.main()
