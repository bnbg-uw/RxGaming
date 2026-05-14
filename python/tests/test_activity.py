from __future__ import annotations

import os
import sys
from pathlib import Path
import tempfile
from types import SimpleNamespace
import types
import unittest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = ROOT / "python"
if str(PYTHON_DIR) not in sys.path:
    sys.path.insert(0, str(PYTHON_DIR))

try:
    import rxgaming_core  # type: ignore  # noqa: F401
except ImportError:
    stub = types.ModuleType("rxgaming_core")

    class StubProjectSettings:
        def __init__(self, *args: object, **kwargs: object) -> None:
            del args, kwargs

    class StubProjectArea:
        def __init__(self, *args: object, **kwargs: object) -> None:
            del args, kwargs

    stub.ProjectSettings = StubProjectSettings
    stub.ProjectArea = StubProjectArea
    stub.load_project_area = lambda path: object()
    stub.save_project_area = lambda project_area, path: Path(path).write_bytes(b"stub")
    sys.modules["rxgaming_core"] = stub

try:
    from PySide6.QtWidgets import QApplication  # type: ignore
except ModuleNotFoundError:
    QApplication = None

if QApplication is not None:
    import activity as activity_module  # noqa: E402
    import persistence as persistence_module  # noqa: E402
    from PySide6.QtGui import QIcon  # type: ignore  # noqa: E402


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestSaveStateActivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])
        activity_module.Activity.test = True

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(activity_module.Activity, "test"):
            delattr(activity_module.Activity, "test")

    def setUp(self) -> None:
        activity_module.Activity._activities = []
        activity_module.Activity._saved_state = {}
        activity_module.Activity._starting = False
        activity_module.Activity._stopping = False
        activity_module.Activity.try_to_save = False

        self.original_get_save_file_name = activity_module.QFileDialog.getSaveFileName
        self.original_get_existing_directory = activity_module.QFileDialog.getExistingDirectory
        self.original_show_message = activity_module.Activity._show_message
        self.original_quit = activity_module.Activity._app.quit
        self.original_write_file = activity_module.SaveStateActivity.write_file

    def tearDown(self) -> None:
        activity_module.QFileDialog.getSaveFileName = self.original_get_save_file_name
        activity_module.QFileDialog.getExistingDirectory = self.original_get_existing_directory
        activity_module.Activity._show_message = self.original_show_message
        activity_module.Activity._app.quit = self.original_quit
        activity_module.SaveStateActivity.write_file = self.original_write_file

        for current in list(activity_module.Activity._activities):
            current.window.close()
        activity_module.Activity._activities = []
        activity_module.Activity._saved_state = {}
        activity_module.Activity._starting = False
        activity_module.Activity._stopping = False

    def test_get_qapplication_applies_resolved_window_icon(self) -> None:
        expected_icon = QIcon(str(ROOT / "icons" / "Icon.ico"))
        original_resolve_icon = activity_module._resolve_application_icon
        original_app = activity_module.QApplication.instance()

        self.assertFalse(expected_icon.isNull())

        try:
            activity_module._resolve_application_icon = lambda: expected_icon
            app = activity_module._get_qapplication()
        finally:
            activity_module._resolve_application_icon = original_resolve_icon

        self.assertEqual(expected_icon.cacheKey(), app.windowIcon().cacheKey())
        if original_app is not None:
            original_app.setWindowIcon(QIcon())

    def test_no_clicked_quits_when_save_prompt_is_last_activity(self) -> None:
        activity = self._make_activity()
        quit_calls: list[str] = []
        activity_module.Activity._stopping = True
        activity_module.Activity._app.quit = lambda: quit_calls.append("quit")

        activity.no_clicked()
        self.app.processEvents()

        self.assertEqual(["quit"], quit_calls)
        self.assertEqual([], activity_module.Activity._activities)
        self.assertFalse(activity_module.Activity._stopping)

    def test_window_close_quits_when_save_prompt_is_last_activity(self) -> None:
        activity = self._make_activity()
        quit_calls: list[str] = []
        activity_module.Activity._stopping = True
        activity_module.Activity._app.quit = lambda: quit_calls.append("quit")

        activity.window.close()
        self.app.processEvents()

        self.assertEqual(["quit"], quit_calls)
        self.assertEqual([], activity_module.Activity._activities)
        self.assertFalse(activity_module.Activity._stopping)

    def test_yes_clicked_stays_open_when_user_cancels_save_dialog(self) -> None:
        activity = self._make_activity()
        messages: list[tuple[str, str]] = []
        activity_module.QFileDialog.getSaveFileName = staticmethod(lambda *args, **kwargs: ("", ""))
        activity_module.Activity._show_message = staticmethod(
            lambda icon, title, text, informative_text=None: messages.append((title, text))
        )

        activity.yes_clicked()
        self.app.processEvents()

        self.assertEqual(1, len(messages))
        self.assertIn(activity, activity_module.Activity._activities)

    def test_yes_clicked_saves_then_quits_when_stopping(self) -> None:
        activity = self._make_activity(
            {
                "ProjectSettingsForm": {
                    "project_name": "Demo",
                }
            }
        )
        quit_calls: list[str] = []
        activity_module.Activity._stopping = True
        activity_module.Activity._app.quit = lambda: quit_calls.append("quit")

        with tempfile.TemporaryDirectory() as tmpdir:
            save_path = Path(tmpdir) / "settings.json"
            activity_module.QFileDialog.getSaveFileName = staticmethod(
                lambda *args, **kwargs: (str(save_path), "JSON files (*.json)")
            )

            activity.yes_clicked()
            self.app.processEvents()

            self.assertEqual(["quit"], quit_calls)
            self.assertTrue(save_path.exists())
            self.assertEqual([], activity_module.Activity._activities)
            self.assertFalse(activity_module.Activity._stopping)

    def test_yes_clicked_uses_project_folder_picker_for_project_snapshot(self) -> None:
        activity = self._make_activity(
            {
                "ProjectSettings": object(),
                "ProjectArea": object(),
                "ProjectSettingsForm": {"project_name": "Demo"},
                "SessionState": {},
            }
        )
        captured_paths: list[str] = []
        activity_module.QFileDialog.getExistingDirectory = staticmethod(lambda *args, **kwargs: "C:/tmp/demo-project")
        activity_module.SaveStateActivity.write_file = staticmethod(lambda path, saved_state: captured_paths.append(str(path)))

        activity.yes_clicked()
        self.app.processEvents()

        self.assertEqual(["C:/tmp/demo-project"], captured_paths)

    def _make_activity(
        self,
        saved_state: dict[str, object] | None = None,
    ) -> activity_module.SaveStateActivity:
        activity = activity_module.SaveStateActivity(None, activity_module.WindowMode.SimultaneousParent)
        activity.on_start(saved_state or {})
        activity.window.show()
        activity_module.Activity._activities.append(activity)
        self.app.processEvents()
        return activity


if __name__ == "__main__":
    unittest.main()


@unittest.skipIf(QApplication is None, "PySide6 is not available in the test runtime")
class TestLoadStateActivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.app = QApplication.instance() or QApplication([])
        activity_module.Activity.test = True

    @classmethod
    def tearDownClass(cls) -> None:
        if hasattr(activity_module.Activity, "test"):
            delattr(activity_module.Activity, "test")

    def setUp(self) -> None:
        activity_module.Activity._activities = []
        activity_module.Activity._saved_state = {}
        self.original_show_message = activity_module.Activity._show_message
        self.original_get_existing_directory = activity_module.QFileDialog.getExistingDirectory
        self.original_get_open_file_name = activity_module.QFileDialog.getOpenFileName
        self.original_read_project_settings_file = persistence_module.read_project_settings_file
        self.original_read_project_snapshot = persistence_module.read_project_snapshot

    def tearDown(self) -> None:
        activity_module.Activity._show_message = self.original_show_message
        activity_module.QFileDialog.getExistingDirectory = self.original_get_existing_directory
        activity_module.QFileDialog.getOpenFileName = self.original_get_open_file_name
        persistence_module.read_project_settings_file = self.original_read_project_settings_file
        persistence_module.read_project_snapshot = self.original_read_project_snapshot
        for current in list(activity_module.Activity._activities):
            current.window.close()
        activity_module.Activity._activities = []
        activity_module.Activity._saved_state = {}

    def test_load_settings_clicked_loads_descriptive_settings_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            settings_path = Path(tmpdir) / "dry-creek-template.json"
            settings_path.write_text("{}", encoding="utf-8")
            activity_module.QFileDialog.getOpenFileName = staticmethod(
                lambda *args, **kwargs: (str(settings_path), "JSON files (*.json)")
            )
            persistence_module.read_project_settings_file = lambda path: SimpleNamespace(
                settings_path=Path(path),
                form_state={"project_name": "Demo"},
            )

            activity = self._make_activity()
            activity.load_settings_clicked()
            self.app.processEvents()

            self.assertEqual({"project_name": "Demo"}, activity_module.Activity._saved_state["ProjectSettingsForm"])
            self.assertEqual(str(settings_path), activity_module.Activity._saved_state["save_file_location"])
            self.assertTrue(activity_module.Activity._saved_state["LoadStateContinue"])

    def test_load_project_clicked_loads_project_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = Path(tmpdir) / "demo-project"
            project_root.mkdir()
            activity_module.QFileDialog.getExistingDirectory = staticmethod(lambda *args, **kwargs: str(project_root))
            persistence_module.read_project_snapshot = lambda path: SimpleNamespace(
                project_root=Path(path),
                project_settings={"name": "Demo"},
                project_area=object(),
                form_state={"project_name": "Demo"},
                session_state={"active_page": 1},
            )

            activity = self._make_activity()
            activity.load_project_clicked()
            self.app.processEvents()

            self.assertEqual(str(project_root), activity_module.Activity._saved_state["ProjectSnapshotPath"])
            self.assertEqual({"active_page": 1}, activity_module.Activity._saved_state["SessionState"])
            self.assertTrue(activity_module.Activity._saved_state["LoadStateContinue"])

    def test_load_settings_clicked_rejects_unsupported_file_type(self) -> None:
        messages: list[tuple[str, str]] = []
        with tempfile.TemporaryDirectory() as tmpdir:
            bad_path = Path(tmpdir) / "legacy.dat"
            bad_path.write_text("", encoding="utf-8")
            activity_module.QFileDialog.getOpenFileName = staticmethod(
                lambda *args, **kwargs: (str(bad_path), "All files (*)")
            )
            activity_module.Activity._show_message = staticmethod(
                lambda icon, title, text, informative_text=None: messages.append((title, text))
            )

            activity = self._make_activity()
            activity.load_settings_clicked()
            self.app.processEvents()

            self.assertEqual([("Unsupported file", "Please choose a settings JSON file or a project folder.")], messages)
            self.assertFalse(activity_module.Activity._saved_state)

    def test_new_clicked_marks_startup_to_continue(self) -> None:
        activity = self._make_activity()

        activity.new_clicked()
        self.app.processEvents()

        self.assertTrue(activity_module.Activity._saved_state["LoadStateContinue"])

    def test_window_close_marks_startup_cancelled(self) -> None:
        activity = self._make_activity()

        activity.window.close()
        self.app.processEvents()

        self.assertFalse(activity_module.Activity._saved_state["LoadStateContinue"])

    def _make_activity(self) -> activity_module.LoadStateActivity:
        activity = activity_module.LoadStateActivity(None, activity_module.WindowMode.SimultaneousParent)
        activity.on_start({})
        activity.window.show()
        activity_module.Activity._activities.append(activity)
        self.app.processEvents()
        return activity
