from __future__ import annotations

import os
import sys
import tempfile
import time
import json
import types
from pathlib import Path
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

from test_support_rxgaming_core import ensure_rxgaming_core_test_module

ensure_rxgaming_core_test_module()

if "gaming_ui" not in sys.modules:
    gaming_ui_stub = types.ModuleType("gaming_ui")

    class StubGamingTabs:
        def __init__(self, *args: object, **kwargs: object) -> None:
            del args, kwargs

    gaming_ui_stub.GamingTabs = StubGamingTabs
    sys.modules["gaming_ui"] = gaming_ui_stub

try:
    from PySide6.QtWidgets import QApplication  # type: ignore
    import PySide6.QtWidgets as qt_widgets  # type: ignore
except ModuleNotFoundError:
    QApplication = None
    qt_widgets = None

if QApplication is not None:
    import activity as activity_module  # noqa: E402
    import gamingactivity as gamingactivity_module  # noqa: E402
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
        self.original_start_project_area_build = projectsettingsactivity_module.start_project_area_build
        self.original_poll_project_area_build = projectsettingsactivity_module.poll_project_area_build
        self.original_finish_project_area_build = projectsettingsactivity_module.finish_project_area_build
        self.original_start_activity = activity_module.Activity.start_activity
        self.original_notify_exception = projectsettingsactivity_module.ProjectSettingsActivity.notify_exception
        self.original_get_save_file_name = qt_widgets.QFileDialog.getSaveFileName

    def tearDown(self) -> None:
        projectsettingsactivity_module.ProjectSettings = self.original_project_settings
        projectsettingsactivity_module.ProjectArea = self.original_project_area
        projectsettingsactivity_module.start_project_area_build = self.original_start_project_area_build
        projectsettingsactivity_module.poll_project_area_build = self.original_poll_project_area_build
        projectsettingsactivity_module.finish_project_area_build = self.original_finish_project_area_build
        activity_module.Activity.start_activity = self.original_start_activity
        projectsettingsactivity_module.ProjectSettingsActivity.notify_exception = self.original_notify_exception
        qt_widgets.QFileDialog.getSaveFileName = self.original_get_save_file_name

    class FakeBuildSnapshot:
        def __init__(
            self,
            *,
            stage: str = "",
            message: str = "",
            completed: int = -1,
            total: int = -1,
            status: str = "running",
            error: str = "",
        ) -> None:
            self.stage = stage
            self.message = message
            self.completed = completed
            self.total = total
            self.status = status
            self.error = error

    class FakeBuildHandle:
        def __init__(
            self,
            snapshots: list["TestProjectSettingsActivity.FakeBuildSnapshot"],
            *,
            result: object | None = None,
            finish_error: Exception | None = None,
        ) -> None:
            self.snapshots = snapshots
            self.result = result
            self.finish_error = finish_error
            self.poll_calls = 0
            self.finish_calls = 0
            self.finalized = False

    def _install_fake_build_api(self, handles: list["TestProjectSettingsActivity.FakeBuildHandle"]) -> list[object]:
        created_handles: list[TestProjectSettingsActivity.FakeBuildHandle] = []

        def fake_start(project_settings: object) -> TestProjectSettingsActivity.FakeBuildHandle:
            del project_settings
            if not handles:
                self.fail("No fake build handles were configured.")
            handle = handles.pop(0)
            created_handles.append(handle)
            return handle

        def fake_poll(handle: TestProjectSettingsActivity.FakeBuildHandle) -> TestProjectSettingsActivity.FakeBuildSnapshot:
            index = min(handle.poll_calls, len(handle.snapshots) - 1)
            handle.poll_calls += 1
            return handle.snapshots[index]

        def fake_finish(handle: TestProjectSettingsActivity.FakeBuildHandle) -> object:
            handle.finish_calls += 1
            if not handle.snapshots:
                raise RuntimeError("No snapshots were configured.")
            terminal = handle.snapshots[min(max(handle.poll_calls - 1, 0), len(handle.snapshots) - 1)]
            if terminal.status == "running":
                raise RuntimeError("Project build is still running.")
            if handle.finalized:
                raise RuntimeError("Project build has already been finalized.")
            handle.finalized = True
            if handle.finish_error is not None:
                raise handle.finish_error
            return handle.result

        projectsettingsactivity_module.start_project_area_build = fake_start
        projectsettingsactivity_module.poll_project_area_build = fake_poll
        projectsettingsactivity_module.finish_project_area_build = fake_finish
        return created_handles

    def test_polling_updates_text_box_before_build_finishes(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        class FakeProjectArea:
            pass

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [
                        self.FakeBuildSnapshot(message="step 1", status="running"),
                        self.FakeBuildSnapshot(message="step 2", status="running"),
                        self.FakeBuildSnapshot(message="finished", status="succeeded"),
                    ],
                    result=FakeProjectArea(),
                )
            ]
        )

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
            self.assertLess(seen_step_one_at - started_at, 0.5)
            self.assertNotIn("step 2", activity.text_output.toPlainText())

            self._wait_until(lambda: bool(launched))
            self.assertTrue(activity.start_button.isEnabled())
            self.assertIn("step 2", activity.text_output.toPlainText())

    def test_duplicate_snapshots_do_not_append_duplicate_lines(self) -> None:
        class FakeProjectArea:
            pass

        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [
                        self.FakeBuildSnapshot(message="Loading lidar", status="running"),
                        self.FakeBuildSnapshot(message="Loading lidar", status="running"),
                        self.FakeBuildSnapshot(message="finished", status="succeeded"),
                    ],
                    result=FakeProjectArea(),
                )
            ]
        )

        launched = []
        activity_module.Activity.start_activity = staticmethod(lambda *args, **kwargs: launched.append((args, kwargs)))

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            activity.start_clicked()
            self._wait_until(lambda: bool(launched))
            self.assertEqual(1, activity.text_output.toPlainText().count("Loading lidar"))

    def test_advancing_completed_with_same_message_updates_visible_progress(self) -> None:
        class FakeProjectArea:
            pass

        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [
                        self.FakeBuildSnapshot(message="Processing units", completed=0, total=10, status="running"),
                        self.FakeBuildSnapshot(message="Processing units", completed=5, total=10, status="running"),
                        self.FakeBuildSnapshot(message="Finished building", completed=10, total=10, status="succeeded"),
                    ],
                    result=FakeProjectArea(),
                )
            ]
        )

        launched = []
        activity_module.Activity.start_activity = staticmethod(lambda *args, **kwargs: launched.append((args, kwargs)))

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            activity.start_clicked()
            self._wait_until(lambda: "[5/10]" in activity.text_output.toPlainText())
            self._wait_until(lambda: bool(launched))
            self.assertIn("Processing units [0/10]", activity.text_output.toPlainText())
            self.assertIn("Processing units [5/10]", activity.text_output.toPlainText())

    def test_success_launches_gaming_activity_with_project_objects(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        class FakeProjectArea:
            pass

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [self.FakeBuildSnapshot(message="finished", status="succeeded")],
                    result=FakeProjectArea(),
                )
            ]
        )

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
            self.assertIsNone(args[2]["ProjectSnapshotPath"])
            self.assertEqual({}, args[2]["SessionState"])
            self.assertEqual("Test Project", args[2]["ProjectSettingsForm"]["project_name"])
            self.assertIs(args[3], projectsettingsactivity_module.WindowMode.Sibling)
            self.assertEqual({}, kwargs)
            self.assertTrue(activity.start_button.isEnabled())

    def test_failure_reenables_start_and_reports_error(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [
                        self.FakeBuildSnapshot(message="before error", status="running"),
                        self.FakeBuildSnapshot(message="before error", status="failed", error="boom"),
                    ],
                    finish_error=RuntimeError("boom"),
                )
            ]
        )

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

    def test_closing_window_stops_polling_and_disposes_handle(self) -> None:
        handle = self.FakeBuildHandle(
            [
                self.FakeBuildSnapshot(message="step 1", status="running"),
                self.FakeBuildSnapshot(message="step 1", status="running"),
            ]
        )
        self._install_fake_build_api([handle])

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            activity.start_clicked()
            self._wait_until(lambda: handle.poll_calls > 0)
            worker = activity.worker
            self.assertIsNotNone(worker)
            activity.window.close()
            self.app.processEvents()
            self.assertFalse(activity._progress_timer.isActive())
            self.assertIsNone(activity.worker)
            assert worker is not None
            self.assertIsNone(worker.handle)

    def test_finalize_raises_while_build_is_running(self) -> None:
        worker = projectsettingsactivity_module.ProjectBuildWorker(
            "project",
            "unit",
            "ref",
            "prop",
            "fia",
            "lidar",
            "name",
            "save",
            4,
        )
        self._install_fake_build_api(
            [self.FakeBuildHandle([self.FakeBuildSnapshot(message="running", status="running")])]
        )

        worker.start()

        with self.assertRaisesRegex(RuntimeError, "still running"):
            worker.finalize()

    def test_finalize_rejects_double_finish(self) -> None:
        class FakeProjectArea:
            pass

        worker = projectsettingsactivity_module.ProjectBuildWorker(
            "project",
            "unit",
            "ref",
            "prop",
            "fia",
            "lidar",
            "name",
            "save",
            4,
        )
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [self.FakeBuildSnapshot(message="done", status="succeeded")],
                    result=FakeProjectArea(),
                )
            ]
        )

        worker.start()
        self.assertIsInstance(worker.finalize(), FakeProjectArea)
        with self.assertRaisesRegex(RuntimeError, "already been finalized"):
            worker.finalize()

    def test_save_settings_as_preserves_descriptive_filename(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            save_path = Path(tmpdir) / "demo-template.json"
            qt_widgets.QFileDialog.getSaveFileName = staticmethod(
                lambda *args, **kwargs: (str(save_path), "RxGaming settings files (*.json)")
            )

            activity.save_as_clicked()

            self.assertEqual(str(save_path), activity.save_file_location)
            payload = json.loads(save_path.read_text(encoding="utf-8"))
            self.assertEqual("rxgaming-settings", payload["format"])
            self.assertEqual("Test Project", payload["form_state"]["project_name"])

    def test_autosave_accepts_new_project_folder_when_parent_exists(self) -> None:
        class FakeProjectSettings:
            def __init__(self, *args: object) -> None:
                self.args = args

        class FakeProjectArea:
            pass

        projectsettingsactivity_module.ProjectSettings = FakeProjectSettings
        projectsettingsactivity_module.ProjectArea = FakeProjectArea
        self._install_fake_build_api(
            [
                self.FakeBuildHandle(
                    [self.FakeBuildSnapshot(message="finished", status="succeeded")],
                    result=FakeProjectArea(),
                )
            ]
        )

        launched = []
        activity_module.Activity.start_activity = staticmethod(lambda *args, **kwargs: launched.append((args, kwargs)))

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            target_folder = Path(tmpdir) / "new-project-folder"
            activity.auto_save_checkbox.setChecked(True)
            activity.auto_save_line_edit.set_text(str(target_folder))
            activity.start_clicked()
            self._wait_until(lambda: len(launched) == 1)

            args, _kwargs = launched[0]
            self.assertEqual(str(target_folder), args[2]["ProjectSnapshotPath"])

    def test_autosave_rejects_missing_parent_folder(self) -> None:
        errors = []
        projectsettingsactivity_module.ProjectSettingsActivity.notify_exception = staticmethod(errors.append)

        with tempfile.TemporaryDirectory() as tmpdir:
            activity = self._make_activity(Path(tmpdir))
            missing_parent = Path(tmpdir) / "missing-parent" / "project-folder"
            activity.auto_save_checkbox.setChecked(True)
            activity.auto_save_line_edit.set_text(str(missing_parent))

            activity.start_clicked()

            self.assertEqual(
                ["The parent folder for the project folder does not exist. Enter a valid folder path before continuing."],
                errors,
            )

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


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestGamingActivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])
        activity_module.Activity.test = True

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(activity_module.Activity, "test"):
            delattr(activity_module.Activity, "test")

    def setUp(self) -> None:
        self.original_project_settings = gamingactivity_module.ProjectSettings
        self.original_project_area = gamingactivity_module.ProjectArea
        self.original_tabs = gamingactivity_module.GamingTabs
        self.original_persistence = gamingactivity_module.ProjectSnapshotSessionPersistence
        self.original_export_georeferenced_raster = gamingactivity_module.export_georeferenced_raster
        self.original_notify_success = gamingactivity_module.GamingActivity._notify_save_success
        self.original_notify_failure = gamingactivity_module.GamingActivity._notify_save_failure
        self.original_question = qt_widgets.QMessageBox.question

    def tearDown(self) -> None:
        gamingactivity_module.ProjectSettings = self.original_project_settings
        gamingactivity_module.ProjectArea = self.original_project_area
        gamingactivity_module.GamingTabs = self.original_tabs
        gamingactivity_module.ProjectSnapshotSessionPersistence = self.original_persistence
        gamingactivity_module.export_georeferenced_raster = self.original_export_georeferenced_raster
        gamingactivity_module.GamingActivity._notify_save_success = self.original_notify_success
        gamingactivity_module.GamingActivity._notify_save_failure = self.original_notify_failure
        qt_widgets.QMessageBox.question = self.original_question
        for widget in self.app.topLevelWidgets():
            widget.close()

    def test_save_project_as_uses_folder_and_followup_save_reuses_it(self) -> None:
        class FakeProjectSettings:
            def __init__(self) -> None:
                self.name = "Demo"
                self.unitPolyPath = "units.shp"
                self.refDataPath = ""
                self.mcsPropPath = "props.csv"
                self.fiaPath = "fia.csv"
                self.lidarPath = "lidar"
                self.unitName = "NAME"
                self.savePath = ""
                self.nThread = 2

        class FakeProjectArea:
            pass

        class FakeSessionState:
            def to_dict(self) -> dict[str, object]:
                return {"active_page": 1}

        class FakeTabs(projectsettingsactivity_module.QTextBrowser):
            def __init__(self, *args: object, **kwargs: object) -> None:
                super().__init__()
                self.session_state = FakeSessionState()

        save_calls: list[str] = []

        class FakePersistence:
            def __init__(self, project_root: str, **kwargs: object) -> None:
                self.project_root = project_root

            def initialize_snapshot(self, state: object) -> None:
                del state

            def save_full_project(self, state: object) -> Path:
                del state
                save_calls.append(self.project_root)
                return Path(self.project_root)

        gamingactivity_module.ProjectSettings = FakeProjectSettings
        gamingactivity_module.ProjectArea = FakeProjectArea
        gamingactivity_module.GamingTabs = FakeTabs
        gamingactivity_module.ProjectSnapshotSessionPersistence = FakePersistence
        gamingactivity_module.GamingActivity._notify_save_success = staticmethod(lambda text: None)
        gamingactivity_module.GamingActivity._notify_save_failure = staticmethod(lambda text: self.fail(text))

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = Path(tmpdir) / "saved-project"

            import PySide6.QtWidgets as qt_widgets

            original_get_existing_directory = qt_widgets.QFileDialog.getExistingDirectory
            qt_widgets.QFileDialog.getExistingDirectory = staticmethod(lambda *args, **kwargs: str(project_root))
            try:
                activity = gamingactivity_module.GamingActivity(None, gamingactivity_module.WindowMode.SimultaneousParent)
                activity.on_start(
                    {
                        "ProjectSettings": FakeProjectSettings(),
                        "ProjectArea": FakeProjectArea(),
                        "ProjectSettingsForm": {"project_name": "Demo"},
                        "SessionState": {},
                    }
                )

                activity.save_project_as()
                activity.save_project()

                self.assertEqual([str(project_root), str(project_root)], save_calls)
                self.assertEqual(str(project_root), activity.project_snapshot_path)
            finally:
                qt_widgets.QFileDialog.getExistingDirectory = original_get_existing_directory

    def test_save_project_as_warns_before_overwriting_existing_project_folder(self) -> None:
        class FakeProjectSettings:
            def __init__(self) -> None:
                self.name = "Demo"
                self.unitPolyPath = "units.shp"
                self.refDataPath = ""
                self.mcsPropPath = "props.csv"
                self.fiaPath = "fia.csv"
                self.lidarPath = "lidar"
                self.unitName = "NAME"
                self.savePath = ""
                self.nThread = 2

        class FakeProjectArea:
            pass

        class FakeSessionState:
            def to_dict(self) -> dict[str, object]:
                return {"active_page": 1}

        class FakeTabs(projectsettingsactivity_module.QTextBrowser):
            def __init__(self, *args: object, **kwargs: object) -> None:
                super().__init__()
                self.session_state = FakeSessionState()

        save_calls: list[str] = []
        question_calls: list[str] = []

        class FakePersistence:
            def __init__(self, project_root: str, **kwargs: object) -> None:
                self.project_root = project_root

            def initialize_snapshot(self, state: object) -> None:
                del state

            def save_full_project(self, state: object) -> Path:
                del state
                save_calls.append(self.project_root)
                return Path(self.project_root)

        gamingactivity_module.ProjectSettings = FakeProjectSettings
        gamingactivity_module.ProjectArea = FakeProjectArea
        gamingactivity_module.GamingTabs = FakeTabs
        gamingactivity_module.ProjectSnapshotSessionPersistence = FakePersistence
        gamingactivity_module.GamingActivity._notify_save_success = staticmethod(lambda text: None)
        gamingactivity_module.GamingActivity._notify_save_failure = staticmethod(lambda text: self.fail(text))

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = Path(tmpdir) / "saved-project"
            project_root.mkdir()
            (project_root / "project.json").write_text("{}", encoding="utf-8")

            original_get_existing_directory = qt_widgets.QFileDialog.getExistingDirectory
            qt_widgets.QFileDialog.getExistingDirectory = staticmethod(lambda *args, **kwargs: str(project_root))
            try:
                activity = gamingactivity_module.GamingActivity(None, gamingactivity_module.WindowMode.SimultaneousParent)
                activity.on_start(
                    {
                        "ProjectSettings": FakeProjectSettings(),
                        "ProjectArea": FakeProjectArea(),
                        "ProjectSettingsForm": {"project_name": "Demo"},
                        "SessionState": {},
                    }
                )

                qt_widgets.QMessageBox.question = staticmethod(
                    lambda *args, **kwargs: question_calls.append(str(project_root)) or qt_widgets.QMessageBox.StandardButton.No
                )
                activity.save_project_as()
                self.assertEqual([str(project_root)], question_calls)
                self.assertEqual([], save_calls)

                qt_widgets.QMessageBox.question = staticmethod(
                    lambda *args, **kwargs: qt_widgets.QMessageBox.StandardButton.Yes
                )
                activity.save_project_as()
                self.assertEqual([str(project_root)], save_calls)
            finally:
                qt_widgets.QFileDialog.getExistingDirectory = original_get_existing_directory

    def test_georeferenced_raster_export_action_is_available(self) -> None:
        class FakeProjectSettings:
            def __init__(self) -> None:
                self.name = "Demo"
                self.unitPolyPath = "units.shp"
                self.refDataPath = ""
                self.mcsPropPath = "props.csv"
                self.fiaPath = "fia.csv"
                self.lidarPath = "lidar"
                self.unitName = "NAME"
                self.savePath = ""
                self.nThread = 2

        class FakeProjectArea:
            pass

        class FakeSessionState:
            def to_dict(self) -> dict[str, object]:
                return {}

        class FakeTabs(projectsettingsactivity_module.QTextBrowser):
            def __init__(self, *args: object, **kwargs: object) -> None:
                super().__init__()
                self.session_state = FakeSessionState()

        triggered: list[object] = []

        gamingactivity_module.ProjectSettings = FakeProjectSettings
        gamingactivity_module.ProjectArea = FakeProjectArea
        gamingactivity_module.GamingTabs = FakeTabs
        gamingactivity_module.export_georeferenced_raster = lambda tabs, window: triggered.append((tabs, window))

        activity = gamingactivity_module.GamingActivity(None, gamingactivity_module.WindowMode.SimultaneousParent)
        activity.on_start(
            {
                "ProjectSettings": FakeProjectSettings(),
                "ProjectArea": FakeProjectArea(),
                "ProjectSettingsForm": {"project_name": "Demo"},
                "SessionState": {},
            }
        )

        self.assertEqual('&Export georeferenced raster image ("*.tif")', activity.export_georeferenced_raster_action.text())
        activity.export_georeferenced_raster()
        self.assertEqual(1, len(triggered))
