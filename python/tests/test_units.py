from __future__ import annotations

import re
import sys
from pathlib import Path
import importlib.util
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
UNITS_PATH = ROOT / "python" / "gaming_ui" / "units.py"
spec = importlib.util.spec_from_file_location("gaming_ui.units", UNITS_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError("Could not load gaming_ui.units")
units = importlib.util.module_from_spec(spec)
spec.loader.exec_module(units)


class TestUnits(unittest.TestCase):
    def test_dbh_round_trip(self) -> None:
        metric_value = 30.48
        display_value = units.dbh_to_display(metric_value, units.UnitSystem.IMPERIAL)
        self.assertAlmostEqual(12.0, display_value)
        self.assertAlmostEqual(metric_value, units.dbh_from_display(display_value, units.UnitSystem.IMPERIAL))

    def test_height_round_trip(self) -> None:
        metric_value = 10.0
        display_value = units.height_to_display(metric_value, units.UnitSystem.IMPERIAL)
        self.assertAlmostEqual(10.0 * units.FEET_PER_METER, display_value)
        self.assertAlmostEqual(metric_value, units.height_from_display(display_value, units.UnitSystem.IMPERIAL))

    def test_tph_round_trip(self) -> None:
        metric_value = 247.105381
        display_value = units.tph_to_display(metric_value, units.UnitSystem.IMPERIAL)
        self.assertAlmostEqual(metric_value * units.TPA_PER_TPH, display_value)
        self.assertAlmostEqual(metric_value, units.tph_from_display(display_value, units.UnitSystem.IMPERIAL))

    def test_ba_round_trip(self) -> None:
        metric_value = 30.0
        display_value = units.ba_to_display(metric_value, units.UnitSystem.IMPERIAL)
        self.assertAlmostEqual(metric_value * units.FT2_PER_ACRE_PER_M2_PER_HA, display_value)
        self.assertAlmostEqual(metric_value, units.ba_from_display(display_value, units.UnitSystem.IMPERIAL))

    def test_cc_round_trip(self) -> None:
        metric_value = 0.32
        display_value = units.to_display("cc", metric_value, units.UnitSystem.IMPERIAL)
        self.assertAlmostEqual(32.0, display_value)
        self.assertAlmostEqual(metric_value, units.from_display("cc", display_value, units.UnitSystem.IMPERIAL))

    def test_cc_array_to_display_uses_percent_scale(self) -> None:
        values = np.asarray([0.1, 0.25, 0.9], dtype=float)
        display_values = units.array_to_display("cc", values, units.UnitSystem.METRIC)
        self.assertEqual([10.0, 25.0, 90.0], list(display_values))

    def test_constants_match_unit_hpp(self) -> None:
        unit_hpp = ROOT / "src" / "lapisgis" / "src" / "Unit.hpp"
        content = unit_hpp.read_text(encoding="utf-8")

        foot = self._extract_constant(content, "internationalFoot")
        inch = self._extract_constant(content, "internationalInch")
        hectare = self._extract_constant(content, "hectare")
        acre = self._extract_constant(content, "acre")

        self.assertAlmostEqual(0.3048, foot)
        self.assertAlmostEqual(0.0254, inch)
        self.assertAlmostEqual(10000.0, hectare)
        self.assertAlmostEqual(4046.8564224, acre)
        self.assertAlmostEqual(units.CENTIMETERS_PER_INCH, inch / 0.01)

    @staticmethod
    def _extract_constant(content: str, name: str) -> float:
        match = re.search(rf"{name}\{{[^,]+,\s*([0-9.]+)\s*\}}", content)
        if match is None:
            raise AssertionError(f"Could not find {name} in Unit.hpp")
        return float(match.group(1))


if __name__ == "__main__":
    unittest.main()
