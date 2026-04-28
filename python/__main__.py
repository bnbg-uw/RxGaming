"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

__main__.py
Application entry point for running RxGaming with a standard PySide6 startup flow.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import traceback
from typing import Any

def _base_path() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[1]


def _first_existing_path(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        "Could not locate any of the required startup resources:\n"
        + "\n".join(str(candidate) for candidate in candidates)
    )


def _iter_native_dll_directories(base_path: Path | None = None) -> list[Path]:
    root = base_path or _base_path()
    candidates: list[Path] = []

    if getattr(sys, "frozen", False):
        meipass = getattr(sys, "_MEIPASS", None)
        if meipass is not None:
            candidates.append(Path(meipass))
        candidates.append(Path(sys.executable).resolve().parent)
    else:
        candidates.extend(
            [
                root / "python",
                root / "build" / "release" / "Release",
                root / "build" / "debug" / "Debug",
            ]
        )

        vcpkg_root = os.environ.get("VCPKG_ROOT")
        triplet = os.environ.get("VCPKG_TARGET_TRIPLET", "x64-windows")
        if vcpkg_root:
            candidates.extend(
                [
                    Path(vcpkg_root) / "installed" / triplet / "bin",
                    Path(vcpkg_root) / "installed" / triplet / "debug" / "bin",
                ]
            )

    unique_existing: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        candidate = candidate.resolve()
        if candidate.exists() and candidate not in seen:
            unique_existing.append(candidate)
            seen.add(candidate)
    return unique_existing


def _bootstrap_native_dll_directories(base_path: Path | None = None) -> list[Any]:
    if os.name != "nt":
        return []
    return [os.add_dll_directory(str(path)) for path in _iter_native_dll_directories(base_path)]


_DLL_DIR_HANDLES = _bootstrap_native_dll_directories()

from activity import Activity, LoadStateActivity
from projectsettingsactivity import ProjectSettingsActivity
from rxgaming_core import set_proj_db_path


def resolve_prop_table_path(base_path: Path | None = None) -> Path:
    root = base_path or _base_path()
    return _first_existing_path(
        root / "resources" / "mcs_prop.csv",
        root / "mcs_prop.csv",
    )


def resolve_fia_path(base_path: Path | None = None) -> Path:
    root = base_path or _base_path()
    return _first_existing_path(
        root / "resources" / "fia",
        root / "fia",
    )


def resolve_proj_data_path(base_path: Path | None = None) -> Path:
    root = base_path or _base_path()
    resources_dir = root / "resources"
    _first_existing_path(resources_dir / "proj.db")
    return resources_dir


def proj_data_directory_arg(base_path: Path | None = None) -> str:
    proj_data_path = resolve_proj_data_path(base_path)
    return f"{proj_data_path}{os.sep}"


class AppContext:
    """Coordinates application startup and activity selection."""

    def __init__(self) -> None:
        self._to_start = ProjectSettingsActivity

    def run(self, **kwargs: Any) -> None:
        Activity.start_activity(LoadStateActivity, saved_state={"onLoad": self.on_load})
        if not Activity._saved_state.get("LoadStateContinue"):
            return
        Activity.try_to_save = True
        Activity.start_activity(self._to_start, **kwargs)

    def on_load(self, saved_state: dict[str, Any]) -> None:
        if "ProjectArea" in saved_state:
            from gamingactivity import GamingActivity

            self._to_start = GamingActivity
            return
        if "ProjectSettingsForm" in saved_state:
            self._to_start = ProjectSettingsActivity
            return
        if "LastActivity" in saved_state:
            self._to_start = saved_state["LastActivity"]


def handle_exception(exc_type: type[BaseException], exc_value: BaseException, exc_traceback: Any) -> None:
    print("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    sys.exit(1)


def main() -> int:
    sys.excepthook = handle_exception

    app_context = AppContext()
    set_proj_db_path(proj_data_directory_arg())
    prop_table_path = resolve_prop_table_path()
    fia_path = resolve_fia_path()
    app_context.run(prop_table_path=str(prop_table_path), fia_path=str(fia_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
