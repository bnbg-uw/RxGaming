from __future__ import annotations

import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from gaming_ui.state import StandViewState, default_stand_view_state_payload, parse_stand_view_state_payload


class TestStandViewState(unittest.TestCase):
    def test_to_dict_and_from_dict_round_trip(self) -> None:
        payload = StandViewState(
            selected_unit_index=2,
            active_page=1,
            raster_mode=3,
            show_treatment=True,
            dbh_cutoff=42.5,
        ).to_dict()

        state = StandViewState.from_dict(payload)

        self.assertEqual(2, state.selected_unit_index)
        self.assertEqual(1, state.active_page)
        self.assertEqual(3, state.raster_mode)
        self.assertTrue(state.show_treatment)
        self.assertEqual(42.5, state.dbh_cutoff)
        self.assertEqual(payload, state.to_dict())

    def test_parse_payload_uses_defaults_for_missing_session_state(self) -> None:
        self.assertEqual(default_stand_view_state_payload(), parse_stand_view_state_payload(None))

    def test_parse_payload_uses_defaults_for_invalid_activity_bag_values(self) -> None:
        payload = parse_stand_view_state_payload(
            {
                "selected_unit_index": "2",
                "active_page": object(),
                "show_treatment": "yes",
                "dbh_cutoff": "42.5",
                "unit_system": "metric",
            }
        )

        self.assertEqual(0, payload["selected_unit_index"])
        self.assertEqual(0, payload["active_page"])
        self.assertEqual(0, payload["raster_mode"])
        self.assertFalse(payload["show_treatment"])
        self.assertEqual(76.2, payload["dbh_cutoff"])
        self.assertEqual("metric", payload["unit_system"])

    def test_parse_payload_accepts_partial_valid_mapping(self) -> None:
        payload = parse_stand_view_state_payload(
            {
                "active_page": 1,
                "show_treatment": True,
            }
        )

        self.assertEqual(0, payload["selected_unit_index"])
        self.assertEqual(1, payload["active_page"])
        self.assertEqual(0, payload["raster_mode"])
        self.assertTrue(payload["show_treatment"])
        self.assertEqual(76.2, payload["dbh_cutoff"])
        self.assertEqual("imperial", payload["unit_system"])

    def test_parse_payload_strict_mode_rejects_invalid_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "selected_unit_index"):
            parse_stand_view_state_payload({"selected_unit_index": "2"}, strict=True)


if __name__ == "__main__":
    unittest.main()
