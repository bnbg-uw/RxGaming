from __future__ import annotations

import os
import sys
import types
from enum import StrEnum
from pathlib import Path
import unittest

import numpy as np
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

gaming_ui_stub = types.ModuleType("gaming_ui")
gaming_ui_stub.GamingTabs = object
gaming_ui_stub.__path__ = []
sys.modules.setdefault("gaming_ui", gaming_ui_stub)

units_stub = types.ModuleType("gaming_ui.units")


class StubUnitSystem(StrEnum):
    METRIC = "metric"
    IMPERIAL = "imperial"


units_stub.UnitSystem = StubUnitSystem
units_stub.array_to_display = lambda _metric_kind, values, _unit_system: values
units_stub.display_name_for = lambda metric_kind, _unit_system: metric_kind
sys.modules.setdefault("gaming_ui.units", units_stub)

try:
    from PySide6.QtWidgets import QApplication  # type: ignore  # noqa: E402
except ModuleNotFoundError:
    QApplication = None

if QApplication is not None:
    import gaming_export as gaming_export_module  # noqa: E402


class FakeUnit:
    def __init__(self) -> None:
        self.name = "Demo Unit"
        self.calls: list[tuple[object, ...]] = []

    def export_rendered_geotiff(self, *args: object) -> None:
        self.calls.append(args)


class FakeTabs:
    def __init__(self, unit: FakeUnit | None, raster_mode: int = 0, show_treatment: bool = False) -> None:
        self._unit = unit
        self._raster_mode = raster_mode
        self._show_treatment = show_treatment

    def current_unit(self) -> FakeUnit:
        if self._unit is None:
            raise IndexError("No units available")
        return self._unit

    def unit_system(self) -> StubUnitSystem:
        return StubUnitSystem.METRIC

    def current_raster_mode(self) -> int:
        return self._raster_mode

    def showing_treatment_view(self) -> bool:
        return self._show_treatment

    def currentIndex(self) -> int:
        return 1


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestGamingExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self) -> None:
        self.original_get_save_file_name = gaming_export_module.QFileDialog.getSaveFileName
        self.original_render = gaming_export_module._render_georeferenced_raster
        self.original_show_warning = gaming_export_module._show_warning

    def tearDown(self) -> None:
        gaming_export_module.QFileDialog.getSaveFileName = self.original_get_save_file_name
        gaming_export_module._render_georeferenced_raster = self.original_render
        gaming_export_module._show_warning = self.original_show_warning

    def test_export_georeferenced_raster_uses_current_variant_from_any_screen(self) -> None:
        unit = FakeUnit()
        tabs = FakeTabs(unit, raster_mode=1, show_treatment=True)
        render_calls: list[tuple[object, object, object]] = []

        gaming_export_module.QFileDialog.getSaveFileName = staticmethod(lambda *args, **kwargs: ("C:/tmp/exported_raster", ""))

        def fake_render(selected_unit: object, option: object, unit_system: object) -> object:
            render_calls.append((selected_unit, option, unit_system))
            return gaming_export_module._RenderedGeoRaster(
                image=(255 * np.ones((12, 10, 4), dtype=np.uint8)),
                map_left_px=1,
                map_top_px=2,
                map_width_px=7,
                map_height_px=8,
            )

        gaming_export_module._render_georeferenced_raster = fake_render

        gaming_export_module.export_georeferenced_raster(tabs, None)

        self.assertEqual(1, len(render_calls))
        self.assertIs(unit, render_calls[0][0])
        self.assertEqual("Treated Basin Map", render_calls[0][1].label)
        self.assertEqual(1, len(unit.calls))
        args = unit.calls[0]
        self.assertEqual("C:/tmp/exported_raster.tif", args[0])
        self.assertEqual((12, 10, 4), args[1].shape)
        self.assertEqual((1, 2, 7, 8), args[2:])

    def test_export_georeferenced_raster_uses_current_clump_rendering(self) -> None:
        unit = FakeUnit()
        tabs = FakeTabs(unit, raster_mode=2, show_treatment=False)

        gaming_export_module.QFileDialog.getSaveFileName = staticmethod(lambda *args, **kwargs: ("C:/tmp/clump.tif", ""))

        gaming_export_module._render_georeferenced_raster = lambda *args, **kwargs: gaming_export_module._RenderedGeoRaster(
            image=np.zeros((6, 5, 4), dtype=np.uint8),
            map_left_px=0,
            map_top_px=0,
            map_width_px=5,
            map_height_px=6,
        )

        gaming_export_module.export_georeferenced_raster(tabs, None)

        self.assertEqual(1, len(unit.calls))
        args = unit.calls[0]
        self.assertEqual("C:/tmp/clump.tif", args[0])
        self.assertEqual((6, 5, 4), args[1].shape)
        self.assertEqual((0, 0, 5, 6), args[2:])

    def test_export_georeferenced_raster_warns_when_no_unit_is_available(self) -> None:
        warnings: list[str] = []
        gaming_export_module._show_warning = lambda parent, text: warnings.append(text)

        gaming_export_module.export_georeferenced_raster(FakeTabs(None), None)

        self.assertEqual(["No unit is available to export."], warnings)

    def test_export_georeferenced_raster_returns_when_save_dialog_is_cancelled(self) -> None:
        unit = FakeUnit()
        tabs = FakeTabs(unit)
        gaming_export_module.QFileDialog.getSaveFileName = staticmethod(lambda *args, **kwargs: ("", ""))

        gaming_export_module.export_georeferenced_raster(tabs, None)

        self.assertEqual([], unit.calls)


if __name__ == "__main__":
    unittest.main()
