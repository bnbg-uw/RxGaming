from __future__ import annotations

import csv
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

try:
    from PySide6.QtWidgets import QApplication  # type: ignore
except ModuleNotFoundError:
    QApplication = None

if QApplication is not None:
    from gaming_ui.landscape import LandscapeReferenceTab  # noqa: E402
    from gaming_ui.state import StandViewState  # noqa: E402
    from gaming_ui.units import TPA_PER_TPH, UnitSystem  # noqa: E402


@dataclass
class FakeStructure:
    tph: float
    ba: float
    mcs: float
    cc: float


class FakeUnit:
    def __init__(
        self,
        name: str,
        current: FakeStructure,
        target: FakeStructure,
        treated: FakeStructure | None,
        simulated: list[FakeStructure],
    ) -> None:
        self.name = name
        self.currentStructure = current
        self.targetStructure = target
        self.treatedStructure = treated
        self._simulated = simulated
        self.simulated_calls = 0

    def get_simulated_structures(self, bb_dbh: float) -> list[FakeStructure]:
        self.simulated_calls += 1
        return list(self._simulated)


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestLandscapeTab(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])

    def test_decision_spaces_are_computed_once_and_reused(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            widget, units = self._make_widget(Path(tmpdir))
            self.assertEqual([1, 1], [unit.simulated_calls for unit in units])

            widget.refresh("second_refresh")
            widget._state.dbh_cutoff = 150.0
            widget.refresh("after_cutoff_change")

            self.assertEqual([1, 1], [unit.simulated_calls for unit in units])

    def test_reference_hulls_use_display_units(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            widget, _ = self._make_widget(Path(tmpdir))
            reference_hull = widget._cached_geometry["reference_hulls"]["ba"][0]
            self.assertAlmostEqual(240.0 * TPA_PER_TPH, float(reference_hull[:, 0].max()), places=4)

    def test_hover_shows_annotation_and_arrows_for_each_metric(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            widget, _ = self._make_widget(Path(tmpdir))
            widget.canvas.draw()

            for metric_index, axis in enumerate(widget._axes):
                event = self._event_for_point(widget._unit_lines[metric_index], axis, 0)
                widget._on_hover(event)

                self.assertTrue(widget._annotations[metric_index].get_visible())
                self.assertEqual("Unit 1", widget._annotations[metric_index].get_text())
                self.assertTrue(widget._target_markers[metric_index].get_visible())
                self.assertTrue(widget._treated_markers[metric_index].get_visible())
                self.assertTrue(widget._current_to_target_arrows[metric_index].get_visible())
                self.assertTrue(widget._target_to_treated_arrows[metric_index].get_visible())
                self.assertTrue(all(patch.get_visible() for patch in widget._unit_hull_patches[metric_index][0]))
                self.assertTrue(all(patch.get_visible() for patch in widget._intersection_patches[metric_index][0]))

    def test_mcs_intersection_geometry_is_cached(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            widget, _ = self._make_widget(Path(tmpdir))
            cached_intersection = widget._cached_geometry["units"][1][1]["intersection"]
            cached_patches = widget._intersection_patches[1][1]

            widget.refresh("cached_geometry_check")

            self.assertIs(cached_intersection, widget._cached_geometry["units"][1][1]["intersection"])
            self.assertIs(cached_patches, widget._intersection_patches[1][1])
            self.assertEqual(len(cached_intersection), len(cached_patches))

    def test_click_selects_the_correct_unit(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            widget, _ = self._make_widget(Path(tmpdir))
            picked: list[int] = []
            widget.set_unit_selected_callback(picked.append)
            widget.canvas.draw()

            event = self._event_for_point(widget._unit_lines[1], widget._axes[1], 1, button=1)
            widget._on_click(event)

            self.assertEqual([1], picked)

    def test_degenerate_decision_space_does_not_crash(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            current = FakeStructure(120.0, 25.0, 4.0, 0.30)
            target = FakeStructure(100.0, 18.0, 3.0, 0.25)
            simulated = [FakeStructure(100.0, 20.0, 2.0, 0.20)] * 4
            reference_path = self._write_reference_csv(Path(tmpdir))
            unit = FakeUnit("Degenerate", current, target, None, simulated)

            widget = LandscapeReferenceTab(
                [unit],
                SimpleNamespace(refDataPath=str(reference_path)),
                StandViewState(dbh_cutoff=76.2, unit_system=UnitSystem.IMPERIAL),
            )
            widget.refresh("degenerate")

            self.assertEqual(1, unit.simulated_calls)
            self.assertEqual(1, len(widget._unit_hull_patches[0]))

    def _make_widget(self, tmpdir: Path) -> tuple[LandscapeReferenceTab, list[FakeUnit]]:
        reference_path = self._write_reference_csv(tmpdir)
        units = [
            FakeUnit(
                "Unit 1",
                FakeStructure(180.0, 28.0, 6.0, 0.32),
                FakeStructure(150.0, 22.0, 5.0, 0.28),
                FakeStructure(140.0, 20.0, 4.0, 0.26),
                [
                    FakeStructure(130.0, 20.0, 3.0, 0.24),
                    FakeStructure(150.0, 24.0, 4.0, 0.27),
                    FakeStructure(180.0, 28.0, 6.0, 0.32),
                    FakeStructure(200.0, 31.0, 8.0, 0.36),
                    FakeStructure(170.0, 26.0, 5.0, 0.30),
                ],
            ),
            FakeUnit(
                "Unit 2",
                FakeStructure(220.0, 35.0, 12.0, 0.44),
                FakeStructure(170.0, 26.0, 7.0, 0.34),
                FakeStructure(165.0, 24.0, 6.0, 0.32),
                [
                    FakeStructure(160.0, 25.0, 5.0, 0.30),
                    FakeStructure(180.0, 28.0, 7.0, 0.35),
                    FakeStructure(200.0, 31.0, 10.0, 0.40),
                    FakeStructure(220.0, 35.0, 12.0, 0.44),
                    FakeStructure(240.0, 38.0, 14.0, 0.48),
                ],
            ),
        ]
        widget = LandscapeReferenceTab(
            units,
            SimpleNamespace(refDataPath=str(reference_path)),
            StandViewState(dbh_cutoff=76.2, unit_system=UnitSystem.IMPERIAL),
        )
        widget.resize(1200, 500)
        widget.show()
        self.app.processEvents()
        widget.refresh("test_setup")
        self.app.processEvents()
        return widget, units

    def _write_reference_csv(self, tmpdir: Path) -> Path:
        path = tmpdir / "reference.csv"
        with path.open("w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(fp, fieldnames=["tph", "ba", "mcs", "cc"])
            writer.writeheader()
            writer.writerows(
                [
                    {"tph": "150", "ba": "20", "mcs": "3", "cc": "0.25"},
                    {"tph": "180", "ba": "25", "mcs": "5", "cc": "0.30"},
                    {"tph": "210", "ba": "30", "mcs": "7", "cc": "0.35"},
                    {"tph": "240", "ba": "34", "mcs": "11", "cc": "0.42"},
                ]
            )
        return path

    @staticmethod
    def _event_for_point(line: object, axis: object, point_index: int, *, button: int | None = None) -> object:
        x_data = line.get_xdata()[point_index]
        y_data = line.get_ydata()[point_index]
        x_pixel, y_pixel = axis.transData.transform((x_data, y_data))
        return SimpleNamespace(
            x=float(x_pixel),
            y=float(y_pixel),
            xdata=float(x_data),
            ydata=float(y_data),
            inaxes=axis,
            button=button,
        )


if __name__ == "__main__":
    unittest.main()
