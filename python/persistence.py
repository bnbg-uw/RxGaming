from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Mapping

from gaming_ui.state import GamingSessionPersistence, StandViewState
from rxgaming_core import ProjectSettings, ProjectArea, load_project_area, save_project_area


SCHEMA_VERSION = 1
PROJECT_MANIFEST_NAME = "project.json"
SETTINGS_FILE_NAME = "settings.json"
SESSION_FILE_NAME = "session.json"
PROJECTAREA_FILE_NAME = "projectarea.h5"


@dataclass
class LoadedProjectSettings:
    settings_path: Path
    form_state: dict[str, Any]


@dataclass
class LoadedProjectSnapshot:
    project_root: Path
    project_settings: ProjectSettings
    project_area: ProjectArea
    session_state: dict[str, Any]
    form_state: dict[str, Any] | None = None


def build_form_state(
    *,
    project_name: str,
    unit_poly_path: str,
    reference_data_path: str,
    lidar_data_path: str,
    unit_name: str,
    threads: int,
    auto_save_enabled: bool,
    auto_save_path: str,
) -> dict[str, Any]:
    return {
        "project_name": project_name,
        "unit_poly_path": unit_poly_path,
        "reference_data_path": reference_data_path,
        "lidar_data_path": lidar_data_path,
        "unit_name": unit_name,
        "threads": int(threads),
        "auto_save_enabled": bool(auto_save_enabled),
        "auto_save_path": auto_save_path,
    }


def serialize_project_settings(project_settings: ProjectSettings) -> dict[str, Any]:
    return {
        "name": project_settings.name,
        "unitPolyPath": project_settings.unitPolyPath,
        "refDataPath": project_settings.refDataPath,
        "mcsPropPath": project_settings.mcsPropPath,
        "fiaPath": project_settings.fiaPath,
        "lidarPath": project_settings.lidarPath,
        "unitName": project_settings.unitName,
        "savePath": project_settings.savePath,
        "nThread": int(project_settings.nThread),
    }


def deserialize_project_settings(payload: Mapping[str, Any]) -> ProjectSettings:
    return ProjectSettings(
        payload["name"],
        payload["unitPolyPath"],
        payload.get("refDataPath", ""),
        payload["mcsPropPath"],
        payload["fiaPath"],
        payload["lidarPath"],
        payload.get("unitName", ""),
        payload.get("savePath", ""),
        int(payload["nThread"]),
    )


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", encoding="utf-8", delete=False, dir=path.parent, suffix=".tmp") as tmp:
        json.dump(payload, tmp, indent=2, sort_keys=True)
        tmp.write("\n")
        temp_path = Path(tmp.name)
    os.replace(temp_path, path)


def _atomic_write_native_snapshot(path: Path, project_area: ProjectArea) -> None:
    if save_project_area is None:
        raise RuntimeError("The installed rxgaming_core module does not support native project snapshots yet.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("wb", delete=False, dir=path.parent, suffix=".tmp") as tmp:
        temp_path = Path(tmp.name)
    save_project_area(project_area, str(temp_path))
    os.replace(temp_path, path)


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fp:
        return json.load(fp)


def _normalize_project_root(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.suffix.lower() == ".json":
        raise ValueError("Project saves now target a project folder, not a JSON file.")
    return candidate


def write_project_settings_file(path: str | Path, form_state: Mapping[str, Any], *, app_version: str) -> Path:
    settings_path = Path(path)
    payload = {
        "format": "rxgaming-settings",
        "schema_version": SCHEMA_VERSION,
        "app_version": app_version,
        "form_state": dict(form_state),
    }
    _atomic_write_json(settings_path, payload)
    return settings_path


def read_project_settings_file(path: str | Path) -> LoadedProjectSettings:
    settings_path = Path(path)
    payload = _read_json(settings_path)
    payload_format = payload.get("format")
    if payload_format == "rxgaming-settings":
        if int(payload.get("schema_version", -1)) != SCHEMA_VERSION:
            raise ValueError("Unsupported RxGaming settings schema version.")
        form_state = payload.get("form_state")
        if not isinstance(form_state, dict):
            raise ValueError("Settings file is missing form_state.")
        return LoadedProjectSettings(settings_path=settings_path, form_state=form_state)

    if payload_format == "rxgaming-project-settings":
        if int(payload.get("schema_version", -1)) != SCHEMA_VERSION:
            raise ValueError("Unsupported RxGaming settings schema version.")
        form_state = payload.get("form_state")
        if isinstance(form_state, dict):
            return LoadedProjectSettings(settings_path=settings_path, form_state=form_state)

        project_settings_payload = payload.get("project_settings")
        if not isinstance(project_settings_payload, dict):
            raise ValueError("Project settings file is missing serialized ProjectSettings.")
        reconstructed_form_state = build_form_state(
            project_name=str(project_settings_payload.get("name", "")),
            unit_poly_path=str(project_settings_payload.get("unitPolyPath", "")),
            reference_data_path=str(project_settings_payload.get("refDataPath", "")),
            lidar_data_path=str(project_settings_payload.get("lidarPath", "")),
            unit_name=str(project_settings_payload.get("unitName", "")),
            threads=int(project_settings_payload.get("nThread", 1)),
            auto_save_enabled=bool(project_settings_payload.get("savePath", "")),
            auto_save_path=str(project_settings_payload.get("savePath", "")),
        )
        return LoadedProjectSettings(settings_path=settings_path, form_state=reconstructed_form_state)

    raise ValueError("Selected file is not an RxGaming settings file.")


def write_project_snapshot(
    project_root: str | Path,
    *,
    app_version: str,
    project_settings: ProjectSettings,
    project_area: ProjectArea,
    session_state: Mapping[str, Any],
    form_state: Mapping[str, Any] | None = None,
) -> Path:
    root = _normalize_project_root(project_root)
    root.mkdir(parents=True, exist_ok=True)

    settings_path = root / SETTINGS_FILE_NAME
    session_path = root / SESSION_FILE_NAME
    projectarea_path = root / PROJECTAREA_FILE_NAME
    manifest_file = root / PROJECT_MANIFEST_NAME

    settings_payload: dict[str, Any] = {
        "format": "rxgaming-project-settings",
        "schema_version": SCHEMA_VERSION,
        "project_settings": serialize_project_settings(project_settings),
    }
    if form_state is not None:
        settings_payload["form_state"] = dict(form_state)

    manifest_payload = {
        "format": "rxgaming-project",
        "schema_version": SCHEMA_VERSION,
        "app_version": app_version,
        "files": {
            "settings": SETTINGS_FILE_NAME,
            "session": SESSION_FILE_NAME,
            "project_area": PROJECTAREA_FILE_NAME,
        },
    }

    _atomic_write_json(settings_path, settings_payload)
    _atomic_write_json(session_path, {
        "format": "rxgaming-session",
        "schema_version": SCHEMA_VERSION,
        "session_state": dict(session_state),
    })
    _atomic_write_native_snapshot(projectarea_path, project_area)
    _atomic_write_json(manifest_file, manifest_payload)
    return root


def read_project_snapshot(path: str | Path) -> LoadedProjectSnapshot:
    if load_project_area is None:
        raise RuntimeError("The installed rxgaming_core module does not support native project snapshots yet.")
    project_root = Path(path)
    if not project_root.exists():
        raise FileNotFoundError(f"Project folder does not exist: {project_root}")
    if not project_root.is_dir():
        raise ValueError("Selected path is not a project folder.")
    manifest_path = project_root / PROJECT_MANIFEST_NAME
    manifest = _read_json(manifest_path)
    if manifest.get("format") != "rxgaming-project":
        raise ValueError("Selected file is not an RxGaming project snapshot.")
    if int(manifest.get("schema_version", -1)) != SCHEMA_VERSION:
        raise ValueError("Unsupported RxGaming project schema version.")

    files = manifest.get("files")
    if not isinstance(files, dict):
        raise ValueError("Project manifest is missing the files section.")

    settings_path = project_root / str(files.get("settings", SETTINGS_FILE_NAME))
    session_path = project_root / str(files.get("session", SESSION_FILE_NAME))
    projectarea_path = project_root / str(files.get("project_area", PROJECTAREA_FILE_NAME))

    if not settings_path.exists():
        raise FileNotFoundError(f"Missing settings file: {settings_path}")
    if not session_path.exists():
        raise FileNotFoundError(f"Missing session file: {session_path}")
    if not projectarea_path.exists():
        raise FileNotFoundError(f"Missing native project snapshot file: {projectarea_path}")

    settings_payload = _read_json(settings_path)
    if int(settings_payload.get("schema_version", -1)) != SCHEMA_VERSION:
        raise ValueError("Unsupported RxGaming settings schema version inside project snapshot.")
    project_settings_payload = settings_payload.get("project_settings")
    if not isinstance(project_settings_payload, dict):
        raise ValueError("Project snapshot is missing serialized ProjectSettings.")

    session_payload = _read_json(session_path)
    session_state = session_payload.get("session_state", {})
    if not isinstance(session_state, dict):
        raise ValueError("Project snapshot session.json is malformed.")

    form_state = settings_payload.get("form_state")
    if form_state is not None and not isinstance(form_state, dict):
        raise ValueError("Project snapshot form_state is malformed.")

    return LoadedProjectSnapshot(
        project_root=project_root,
        project_settings=deserialize_project_settings(project_settings_payload),
        project_area=load_project_area(str(projectarea_path)),
        session_state=session_state,
        form_state=form_state,
    )


class ProjectSnapshotSessionPersistence(GamingSessionPersistence):
    def __init__(
        self,
        project_root: str | Path,
        *,
        app_version: str,
        project_settings: ProjectSettings,
        project_area: ProjectArea,
        form_state: Mapping[str, Any] | None = None,
    ) -> None:
        self.project_root = _normalize_project_root(project_root)
        self.app_version = app_version
        self.project_settings = project_settings
        self.project_area = project_area
        self.form_state = dict(form_state) if form_state is not None else None

    def load_initial_state(self, saved_state: dict[str, object]) -> StandViewState:
        return StandViewState.from_saved_state(saved_state)

    def initialize_snapshot(self, state: StandViewState) -> None:
        write_project_snapshot(
            self.project_root,
            app_version=self.app_version,
            project_settings=self.project_settings,
            project_area=self.project_area,
            session_state=state.to_dict(),
            form_state=self.form_state,
        )

    def save_session(self, state: StandViewState, reason: str = "session_updated") -> None:
        session_path = self.project_root / SESSION_FILE_NAME
        _atomic_write_json(session_path, {
            "format": "rxgaming-session",
            "schema_version": SCHEMA_VERSION,
            "session_state": state.to_dict(),
        })
        if reason in {"targets_changed", "show_treatment_toggled", "manual_save"}:
            self.initialize_snapshot(state)

    def save_full_project(self, state: StandViewState) -> Path:
        return write_project_snapshot(
            self.project_root,
            app_version=self.app_version,
            project_settings=self.project_settings,
            project_area=self.project_area,
            session_state=state.to_dict(),
            form_state=self.form_state,
        )
