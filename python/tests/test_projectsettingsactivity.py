from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path
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
    import activity as activity_module  # noqa: E402
    import projectsettingsactivity as projectsettingsactivity_module  # noqa: E402


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestProjectSettingsActivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])
        activity_module.Activity.test = True

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(activity_module.Activity, "test"):
            delattr(activity_module.Activity, "test")

    def setUp(self) -> None:
        self.original_project_settings = projectsettingsactivity_module.ProjectSettings
        self.original_project_area = projectsettingsactivity_module.ProjectArea
        self.original_build_project_area = projectsettingsactivity_module.build_project_area_with_progress
        self.original_start_activity = activity_module.Activity.start_activity
        self.original_notify_exception = projectsettingsactivity_module.ProjectSettingsActivity.notify_exception

    def tearDown(self) -> None:
        projectsettingsactivity_module.ProjectSettings = self.original_project_settings
        projectsettingsactivity_module.ProjectArea = self.original_project_area
        projectsettingsactivity_module.build_project_area_with_progress = self.original_build_project_area
        activity_module.Activity.start_activity = self.original_start_activity
        projectsettingsactivity_module.ProjectSettingsActivity.notify_exception = self.original_notify_exception

    def test_streaming_updates_text_box_before_worker_finishes(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        class FakeProjectArea:
            pass

        class FakeProgressEvent:
            def __init__(self, message: str) -> None:
                self.message = message

        def fake_build_project_area(project_settings: object, callback: object) -> FakeProjectArea:
            del project_settings
            callback(FakeProgressEvent("step 1"))
            time.sleep(0.2)
            callback(FakeProgressEvent("step 2"))
            time.sleep(0.2)
            return FakeProjectArea()

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        projectsettingsactivity_module.build_project_area_with_progress = fake_build_project_area

        launched = []
        activity_module.Activity.start_activity = staticmethod(lambda *args, **kwargs: launched.append((args, kwargs)))

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            started_at = time.monotonic()
            activity.start_clicked()

            self.assertFalse(activity.start_button.isEnabled())
            self._wait_until(lambda: "step 1" in activity.text_output.toPlainText())
            seen_step_one_at = time.monotonic()

            self.assertIn("Processing started:", activity.text_output.toPlainText())
            self.assertIn("step 1", activity.text_output.toPlainText())
            self.assertLess(seen_step_one_at - started_at, 0.35)
            self.assertNotIn("step 2", activity.text_output.toPlainText())

            self._wait_until(lambda: bool(launched))
            self.assertTrue(activity.start_button.isEnabled())
            self.assertIn("step 2", activity.text_output.toPlainText())

    def test_success_launches_gaming_activity_with_project_objects(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        class FakeProjectArea:
            pass

        class FakeProgressEvent:
            def __init__(self, message: str) -> None:
                self.message = message

        def fake_build_project_area(project_settings: object, callback: object) -> FakeProjectArea:
            callback(FakeProgressEvent("finished"))
            self.assertIsInstance(project_settings, FakeProjectSettings)
            return FakeProjectArea()

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        projectsettingsactivity_module.build_project_area_with_progress = fake_build_project_area

        launched = []
        activity_module.Activity.start_activity = staticmethod(lambda *args, **kwargs: launched.append((args, kwargs)))

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            activity.start_clicked()
            self._wait_until(lambda: len(launched) == 1)

            args, kwargs = launched[0]
            self.assertIs(args[0], projectsettingsactivity_module.GamingActivity)
            self.assertIsNone(args[1])
            self.assertIsInstance(args[2]["ProjectSettings"], FakeProjectSettings)
            self.assertIsInstance(args[2]["ProjectArea"], FakeProjectArea)
            self.assertIsNone(args[2]["Autosave_path"])
            self.assertIs(args[3], projectsettingsactivity_module.WindowMode.Sibling)
            self.assertEqual({}, kwargs)
            self.assertTrue(activity.start_button.isEnabled())

    def test_failure_reenables_start_and_reports_error(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        class FakeProgressEvent:
            def __init__(self, message: str) -> None:
                self.message = message

        def fake_build_project_area(project_settings: object, callback: object) -> object:
            del project_settings
            callback(FakeProgressEvent("before error"))
            time.sleep(0.1)
            raise RuntimeError("boom")

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        projectsettingsactivity_module.build_project_area_with_progress = fake_build_project_area

        launched = []
        activity_module.Activity.start_activity = staticmethod(lambda *args, **kwargs: launched.append((args, kwargs)))
        errors = []
        projectsettingsactivity_module.ProjectSettingsActivity.notify_exception = staticmethod(errors.append)

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            activity.start_clicked()

            self._wait_until(lambda: bool(errors))

            self.assertEqual([], launched)
            self.assertTrue(activity.start_button.isEnabled())
            self.assertIn("before error", activity.text_output.toPlainText())
            self.assertIn("boom", errors[0])
            self.assertIn("Debugging Traceback", errors[0])

    def _make_activity(self, tmpdir: Path) -> projectsettingsactivity_module.ProjectSettingsActivity:
        unit_path = tmpdir / "units.shp"
        unit_path.write_text("", encoding="utf-8")
        lidar_path = tmpdir / "lidar"
        lidar_path.mkdir()

        activity = projectsettingsactivity_module.ProjectSettingsActivity(
            None,
            projectsettingsactivity_module.WindowMode.SimultaneousParent,
        )
        activity.on_start({}, prop_table_path="props.csv", fia_path="fia.csv")
        activity.prj_name_edit.setText("Test Project")
        activity.unit_poly_path_edit.set_text(str(unit_path))
        activity.lidar_data_path_edit.set_text(str(lidar_path))
        activity.reference_data_path_edit.set_text("")
        activity.unit_name_edit.setText("UnitName")
        activity.threads_edit.setValue(2)
        activity.window.show()
        self.app.processEvents()
        return activity

    def _wait_until(self, predicate: object, timeout: float = 2.0) -> None:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            self.app.processEvents()
            if predicate():
                return
            time.sleep(0.01)
        self.fail("Timed out waiting for condition")


if __name__ == "__main__":
    unittest.main()
