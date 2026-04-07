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

from pathlib import Path
import sys
import traceback
from typing import Any

from activity import Activity, LoadStateActivity
from projectsettingsactivity import ProjectSettingsActivity


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


def resolve_prop_table_path(base_path: Path | None = None) -> Path:
    root = base_path or _base_path()
    return _first_existing_path(
        root / "resources" / "mcs_prop.csv",
        root / "mcs_prop.csv",
    )


class AppContext:
    """Coordinates application startup and activity selection."""

    def __init__(self) -> None:
        self._to_start = ProjectSettingsActivity

    def run(self, **kwargs: Any) -> None:
        Activity.start_activity(LoadStateActivity, saved_state={"onLoad": self.on_load})
        Activity.try_to_save = True
        Activity.start_activity(self._to_start, **kwargs)

    def on_load(self, saved_state: dict[str, Any]) -> None:
        if "LastActivity" in saved_state:
            self._to_start = saved_state["LastActivity"]


def handle_exception(exc_type: type[BaseException], exc_value: BaseException, exc_traceback: Any) -> None:
    print("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    sys.exit(1)


def main() -> int:
    sys.excepthook = handle_exception

    app_context = AppContext()
    prop_table_path = resolve_prop_table_path()
    app_context.run(prop_table_path=str(prop_table_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
