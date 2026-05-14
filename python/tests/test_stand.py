from __future__ import annotations

import os
from pathlib import Path
import sys
import types
import unittest
import numpy as np

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from test_support_rxgaming_core import ensure_rxgaming_core_test_module

ensure_rxgaming_core_test_module()

from gaming_ui.stand import StandViewCoordinator  # noqa: E402
from gaming_ui.units import UnitSystem  # noqa: E402


class TestStandHelpers(unittest.TestCase):
    def test_tree_dbh_distribution_uses_display_units_and_filters_invalid_values(self) -> None:
        coordinator = StandViewCoordinator.__new__(StandViewCoordinator)
        coordinator.unit_system = UnitSystem.IMPERIAL

        points = np.asarray(
            [
                [0.0, 0.0, 0.0, 0.0, 30.48],
                [0.0, 0.0, 0.0, 0.0, 60.96],
                [0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, np.nan],
            ],
            dtype=float,
        )

        distribution = coordinator._tree_dbh_distribution(points)

        self.assertEqual(2, distribution.size)
        self.assertAlmostEqual(12.0, float(distribution[0]), places=6)
        self.assertAlmostEqual(24.0, float(distribution[1]), places=6)

    def test_tree_dbh_distribution_returns_metric_values_when_metric_is_selected(self) -> None:
        coordinator = StandViewCoordinator.__new__(StandViewCoordinator)
        coordinator.unit_system = UnitSystem.METRIC

        points = np.asarray([[0.0, 0.0, 0.0, 0.0, 25.0]], dtype=float)

        distribution = coordinator._tree_dbh_distribution(points)

        self.assertEqual([25.0], list(distribution))

    def test_cut_ba_per_area_returns_ba_per_acre_in_imperial(self) -> None:
        coordinator = StandViewCoordinator.__new__(StandViewCoordinator)
        coordinator.unit_system = UnitSystem.IMPERIAL

        cut_ba = np.asarray([10.0, 15.0], dtype=float)
        one_acre_in_ha = 4046.8564224 / 10000.0

        ba_per_area = coordinator._cut_ba_per_area(cut_ba, one_acre_in_ha)

        self.assertAlmostEqual(25.0, ba_per_area, places=6)

    def test_format_structure_summary_omits_stand_area(self) -> None:
        structure = types.SimpleNamespace(tph=247.105381, ba=30.0, mcs=4.0, cc=0.5)

        summary = StandViewCoordinator._format_structure_summary(structure, UnitSystem.IMPERIAL)

        self.assertNotIn("Stand Area:", summary)
        self.assertIn("BA: 130.68", summary)


if __name__ == "__main__":
    unittest.main()
