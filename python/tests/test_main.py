from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from test_support_rxgaming_core import ensure_rxgaming_core_test_module

ensure_rxgaming_core_test_module()

try:
    from PySide6.QtWidgets import QApplication  # type: ignore
except ModuleNotFoundError:
    QApplication = None

if QApplication is not None:
    import activity as activity_module  # noqa: E402
    import loadstateactivity as loadstateactivity_module  # noqa: E402

    MAIN_SPEC = importlib.util.spec_from_file_location("rxgaming_main", PYTHON_DIR / "__main__.py")
    assert MAIN_SPEC is not None
    assert MAIN_SPEC.loader is not None
    main_module = importlib.util.module_from_spec(MAIN_SPEC)
    MAIN_SPEC.loader.exec_module(main_module)


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestAppContext(unittest.TestCase):
    def setUp(self) -> None:
        self.original_start_activity = activity_module.Activity.start_activity
        activity_module.Activity._saved_state = {}
        activity_module.Activity.try_to_save = False

    def tearDown(self) -> None:
        activity_module.Activity.start_activity = self.original_start_activity
        activity_module.Activity._saved_state = {}
        activity_module.Activity.try_to_save = False

    def test_run_stops_after_load_state_window_close(self) -> None:
        start_calls: list[object] = []

        def fake_start_activity(activity_class, *args, **kwargs):
            start_calls.append(activity_class)
            if activity_class is loadstateactivity_module.LoadStateActivity:
                activity_module.Activity._saved_state = {"LoadStateContinue": False}
            return None

        activity_module.Activity.start_activity = staticmethod(fake_start_activity)

        context = main_module.AppContext()
        context.run(prop_table_path="prop.csv", fia_path="fia")

        self.assertEqual([loadstateactivity_module.LoadStateActivity], start_calls)
        self.assertFalse(activity_module.Activity.try_to_save)

    def test_run_starts_project_settings_after_new_project_choice(self) -> None:
        start_calls: list[object] = []

        def fake_start_activity(activity_class, *args, **kwargs):
            start_calls.append(activity_class)
            if activity_class is loadstateactivity_module.LoadStateActivity:
                activity_module.Activity._saved_state = {"LoadStateContinue": True}
            return None

        activity_module.Activity.start_activity = staticmethod(fake_start_activity)

        context = main_module.AppContext()
        context.run(prop_table_path="prop.csv", fia_path="fia")

        self.assertEqual(
            [loadstateactivity_module.LoadStateActivity, main_module.ProjectSettingsActivity],
            start_calls,
        )
        self.assertTrue(activity_module.Activity.try_to_save)
