from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace
import unittest

from test_support_rxgaming_core import ensure_rxgaming_core_test_module

ensure_rxgaming_core_test_module()

import persistence  # noqa: E402
from gaming_ui.state import StandViewState  # noqa: E402


class FakeProjectSettings:
    def __init__(
        self,
        name: str,
        unitPolyPath: str,
        refDataPath: str,
        mcsPropPath: str,
        fiaPath: str,
        lidarPath: str,
        unitName: str,
        savePath: str,
        nThread: int,
    ) -> None:
        self.name = name
        self.unitPolyPath = unitPolyPath
        self.refDataPath = refDataPath
        self.mcsPropPath = mcsPropPath
        self.fiaPath = fiaPath
        self.lidarPath = lidarPath
        self.unitName = unitName
        self.savePath = savePath
        self.nThread = nThread


class TestPersistence(unittest.TestCase):
    def setUp(self) -> None:
        self.original_project_settings = persistence.ProjectSettings
        self.original_save_project_area = persistence.save_project_area
        self.original_load_project_area = persistence.load_project_area
        persistence.ProjectSettings = FakeProjectSettings

    def tearDown(self) -> None:
        persistence.ProjectSettings = self.original_project_settings
        persistence.save_project_area = self.original_save_project_area
        persistence.load_project_area = self.original_load_project_area

    def test_settings_file_round_trip(self) -> None:
        form_state = persistence.build_form_state(
            project_name="Demo",
            unit_poly_path="units.shp",
            reference_data_path="reference.csv",
            lidar_data_path="lidar",
            unit_name="NAME",
            threads=4,
            auto_save_enabled=True,
            auto_save_path="C:/tmp/project-folder",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "demo-template.json"
            persistence.write_project_settings_file(path, form_state, app_version="1.2.3")

            loaded = persistence.read_project_settings_file(path)
            self.assertEqual(path, loaded.settings_path)
            self.assertEqual(form_state, loaded.form_state)

    def test_project_snapshot_round_trip_uses_native_payload(self) -> None:
        saved_native_payloads: list[tuple[object, Path]] = []
        fake_project_area = object()

        def fake_save_project_area(project_area: object, path: str) -> None:
            saved_native_payloads.append((project_area, Path(path)))
            Path(path).write_bytes(b"native-snapshot")

        def fake_load_project_area(path: str) -> object:
            self.assertEqual(b"native-snapshot", Path(path).read_bytes())
            return fake_project_area

        persistence.save_project_area = fake_save_project_area
        persistence.load_project_area = fake_load_project_area

        project_settings = FakeProjectSettings(
            "Demo",
            "units.shp",
            "reference.csv",
            "props.csv",
            "fia",
            "missing-lidar-folder",
            "NAME",
            "project-folder",
            8,
        )
        form_state = {"project_name": "Demo", "threads": 8}
        session_state = StandViewState(selected_unit_index=2, active_page=1, raster_mode=3, show_treatment=True).to_dict()

        with tempfile.TemporaryDirectory() as tmpdir:
            project_root = Path(tmpdir) / "project-folder"
            written_root = persistence.write_project_snapshot(
                project_root,
                app_version="2.0.0",
                project_settings=project_settings,
                project_area=fake_project_area,
                session_state=session_state,
                form_state=form_state,
            )

            self.assertEqual(project_root, written_root)
            self.assertEqual(1, len(saved_native_payloads))
            self.assertIs(fake_project_area, saved_native_payloads[0][0])
            self.assertTrue((project_root / persistence.PROJECTAREA_FILE_NAME).exists())

            loaded = persistence.read_project_snapshot(project_root)
            self.assertEqual(project_root, loaded.project_root)
            self.assertIs(fake_project_area, loaded.project_area)
            self.assertEqual(session_state, loaded.session_state)
            self.assertEqual(form_state, loaded.form_state)
            self.assertEqual("missing-lidar-folder", loaded.project_settings.lidarPath)

            manifest_payload = json.loads((project_root / persistence.PROJECT_MANIFEST_NAME).read_text(encoding="utf-8"))
            self.assertEqual("rxgaming-project", manifest_payload["format"])

    def test_project_snapshot_read_requires_project_folder(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            invalid_path = Path(tmpdir) / "project.json"
            invalid_path.write_text("{}", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "project folder"):
                persistence.read_project_snapshot(invalid_path)

    def test_project_settings_file_inside_project_folder_loads_as_settings(self) -> None:
        project_settings = FakeProjectSettings(
            "Demo",
            "units.shp",
            "reference.csv",
            "props.csv",
            "fia",
            "lidar-folder",
            "NAME",
            "project-folder",
            8,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            settings_path = Path(tmpdir) / "settings.json"
            payload = {
                "format": "rxgaming-project-settings",
                "schema_version": persistence.SCHEMA_VERSION,
                "project_settings": persistence.serialize_project_settings(project_settings),
            }
            settings_path.write_text(json.dumps(payload), encoding="utf-8")

            loaded = persistence.read_project_settings_file(settings_path)

            self.assertEqual("Demo", loaded.form_state["project_name"])
            self.assertEqual("units.shp", loaded.form_state["unit_poly_path"])
            self.assertEqual("project-folder", loaded.form_state["auto_save_path"])
            self.assertTrue(loaded.form_state["auto_save_enabled"])

    def test_snapshot_session_persistence_updates_session_file_without_rewriting_snapshot(self) -> None:
        saved_snapshots: list[dict[str, object]] = []

        def fake_write_project_snapshot(
            project_root: str | Path,
            *,
            app_version: str,
            project_settings: object,
            project_area: object,
            session_state: dict[str, object],
            form_state: dict[str, object] | None = None,
        ) -> Path:
            saved_snapshots.append(
                {
                    "project_root": str(project_root),
                    "app_version": app_version,
                    "project_settings": project_settings,
                    "project_area": project_area,
                    "session_state": dict(session_state),
                    "form_state": form_state,
                }
            )
            root = Path(project_root)
            root.mkdir(parents=True, exist_ok=True)
            (root / persistence.PROJECT_MANIFEST_NAME).write_text("{}", encoding="utf-8")
            (root / persistence.SESSION_FILE_NAME).write_text("{}", encoding="utf-8")
            return root

        original_write_project_snapshot = persistence.write_project_snapshot
        persistence.write_project_snapshot = fake_write_project_snapshot
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                project_root = Path(tmpdir) / "project-folder"
                session_persistence = persistence.ProjectSnapshotSessionPersistence(
                    project_root,
                    app_version="9.9.9",
                    project_settings=SimpleNamespace(),
                    project_area=SimpleNamespace(),
                    form_state={"project_name": "Autosave"},
                )
                state = StandViewState(selected_unit_index=1, active_page=2, raster_mode=1, show_treatment=False)

                session_persistence.initialize_snapshot(state)
                session_persistence.save_session(state, "page_changed")

                session_payload = json.loads((project_root / persistence.SESSION_FILE_NAME).read_text(encoding="utf-8"))
                self.assertEqual(state.to_dict(), session_payload["session_state"])
                self.assertEqual(1, len(saved_snapshots))

                session_persistence.save_session(state, "targets_changed")
                session_persistence.save_session(state, "show_treatment_toggled")
                self.assertEqual(1, len(saved_snapshots))
        finally:
            persistence.write_project_snapshot = original_write_project_snapshot


if __name__ == "__main__":
    unittest.main()
