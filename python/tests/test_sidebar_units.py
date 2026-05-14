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

from test_support_rxgaming_core import ensure_rxgaming_core_test_module

ensure_rxgaming_core_test_module()

try:
    from PySide6.QtCore import Qt  # type: ignore
    from PySide6.QtWidgets import QApplication  # type: ignore
except ModuleNotFoundError:
    Qt = None
    QApplication = None

if QApplication is not None:
    from gaming_ui.sidebar import StructureInfo, UnitSidebar, UnitListModel  # noqa: E402
    from gaming_ui.stand import StandViewCoordinator  # noqa: E402
    from gaming_ui.state import StandViewState  # noqa: E402
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
    def __init__(self, name: str = "Unit 1") -> None:
        self.name = name
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

    def test_sidebar_search_filters_units_by_case_insensitive_partial_name(self) -> None:
        sidebar = UnitSidebar(
            [FakeUnit("North Ridge"), FakeUnit("south slope"), FakeUnit("East Basin")],
            units.UnitSystem.IMPERIAL,
        )

        sidebar.unit_search.setText("SoUtH")
        self.app.processEvents()

        self.assertEqual(1, sidebar.filtered_model.rowCount())
        self.assertEqual("south slope", sidebar.filtered_model.index(0, 0).data())

    def test_sidebar_search_clear_restores_full_list(self) -> None:
        sidebar = UnitSidebar(
            [FakeUnit("North Ridge"), FakeUnit("South Slope"), FakeUnit("East Basin")],
            units.UnitSystem.IMPERIAL,
        )

        sidebar.unit_search.setText("ridge")
        self.app.processEvents()
        sidebar.unit_search.clear()
        self.app.processEvents()

        self.assertEqual(3, sidebar.filtered_model.rowCount())

    def test_refresh_preserves_search_text_and_filter(self) -> None:
        sidebar = UnitSidebar(
            [FakeUnit("North Ridge"), FakeUnit("South Slope"), FakeUnit("East Basin")],
            units.UnitSystem.IMPERIAL,
        )

        sidebar.unit_search.setText("east")
        self.app.processEvents()
        sidebar.model.refresh()
        self.app.processEvents()

        self.assertEqual("east", sidebar.unit_search.text())
        self.assertEqual(1, sidebar.filtered_model.rowCount())
        self.assertEqual("East Basin", sidebar.filtered_model.index(0, 0).data())

    def test_filtered_selection_updates_selected_source_index(self) -> None:
        coordinator = self._make_selection_coordinator(
            [FakeUnit("North Ridge"), FakeUnit("South Slope"), FakeUnit("East Basin")]
        )

        coordinator.sidebar.unit_search.setText("east")
        self.app.processEvents()
        coordinator.sidebar.unit_list_view.setCurrentIndex(coordinator.sidebar.filtered_model.index(0, 0))
        self.app.processEvents()

        self.assertEqual(2, coordinator.state.selected_unit_index)

    def test_hidden_selection_is_preserved_when_filter_hides_selected_unit(self) -> None:
        coordinator = self._make_selection_coordinator(
            [FakeUnit("North Ridge"), FakeUnit("South Slope"), FakeUnit("East Basin")]
        )
        coordinator.state.selected_unit_index = 1
        coordinator.select_unit(1)
        self.app.processEvents()

        coordinator.sidebar.unit_search.setText("east")
        self.app.processEvents()

        self.assertEqual(1, coordinator.state.selected_unit_index)
        self.assertFalse(coordinator.sidebar.unit_list_view.currentIndex().isValid())

    def test_selected_unit_is_reselected_when_filter_makes_it_visible_again(self) -> None:
        coordinator = self._make_selection_coordinator(
            [FakeUnit("North Ridge"), FakeUnit("South Slope"), FakeUnit("East Basin")]
        )
        coordinator.state.selected_unit_index = 1
        coordinator.select_unit(1)
        self.app.processEvents()

        coordinator.sidebar.unit_search.setText("east")
        self.app.processEvents()
        coordinator.sidebar.unit_search.clear()
        self.app.processEvents()

        current = coordinator.sidebar.unit_list_view.currentIndex()
        self.assertTrue(current.isValid())
        self.assertEqual("South Slope", current.data())

    def _make_selection_coordinator(self, rx_units: list[FakeUnit]) -> StandViewCoordinator:
        coordinator = StandViewCoordinator.__new__(StandViewCoordinator)
        coordinator.rx_units = rx_units
        coordinator.state = StandViewState(selected_unit_index=0, unit_system=units.UnitSystem.IMPERIAL)
        coordinator.unit_system = coordinator.state.unit_system
        coordinator.sidebar = UnitSidebar(rx_units, coordinator.unit_system)
        coordinator._syncing_sidebar_selection = False
        coordinator._mark_all_pages_dirty = lambda: None
        coordinator.refresh_all = lambda trigger="refresh_all": None
        coordinator._notify_state_changed = lambda reason: None
        coordinator._notify_landscape_invalidated = lambda: None
        coordinator.sidebar.unit_list_view.selectionModel().currentChanged.connect(coordinator._on_unit_changed)
        coordinator.sidebar.unit_search.textChanged.disconnect(coordinator.sidebar.set_unit_filter_text)
        coordinator.sidebar.unit_search.textChanged.connect(coordinator._on_unit_filter_changed)
        return coordinator


if __name__ == "__main__":
    unittest.main()
