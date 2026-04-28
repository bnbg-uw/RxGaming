from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from PySide6.QtWidgets import QFileDialog, QHBoxLayout, QLabel, QMessageBox, QPushButton, QVBoxLayout

from activity import Activity, SavedState

class SaveStateActivity(Activity):
    """Prompts user to store the state instance dictionary."""

    def on_start(self, saved_state: SavedState, **kwargs: Any) -> None:
        del kwargs
        self.label = QLabel("Would you like to save your work?")
        self.yes_button = QPushButton("Yes")
        self.no_button = QPushButton("No")

        self.vlayout = QVBoxLayout()
        self.hlayout = QHBoxLayout()
        self.hlayout.addWidget(self.yes_button)
        self.hlayout.addWidget(self.no_button)
        self.vlayout.addWidget(self.label)
        self.vlayout.addLayout(self.hlayout)
        self.window.setLayout(self.vlayout)
        self.window.setWindowTitle("Save your work")

        self.yes_button.clicked.connect(self.yes_clicked)
        self.no_button.clicked.connect(self.no_clicked)

        self.saved_state = saved_state

    def save(self) -> SavedState:
        return {}

    def yes_clicked(self, checked: bool = False) -> None:
        del checked
        if self.prompt_and_save(self.saved_state):
            self.stop()

    def prompt_and_save(self, saved_state: Mapping[str, Any]) -> bool:
        has_project_snapshot = "ProjectArea" in saved_state and "ProjectSettings" in saved_state
        if has_project_snapshot:
            selected_path = QFileDialog.getExistingDirectory(
                None,
                "Save project folder as...",
                str(self._default_project_folder(saved_state)),
            )
            if selected_path == "":
                Activity._show_message(
                    QMessageBox.Icon.Warning,
                    "Choose a folder",
                    "Please select a project folder to save your work",
                )
                return False
            type(self).write_file(selected_path, saved_state)
            return True

        file_path = QFileDialog.getSaveFileName(
            None,
            "Save settings as...",
            str(self._default_settings_file(saved_state)),
            "RxGaming settings files (*.json)",
        )[0]
        if file_path == "":
            Activity._show_message(
                QMessageBox.Icon.Warning,
                "Choose a file",
                "Please select a settings JSON file to save your work",
            )
            return False
        type(self).write_file(file_path, saved_state)
        return True

    @staticmethod
    def _default_project_folder(saved_state: Mapping[str, Any]) -> Path:
        project_root = saved_state.get("ProjectSnapshotPath")
        if isinstance(project_root, str) and project_root:
            return Path(project_root)

        project_settings = saved_state.get("ProjectSettings")
        save_path = getattr(project_settings, "savePath", "")
        if isinstance(save_path, str) and save_path:
            return Path(save_path)
        return Path.cwd()

    @staticmethod
    def _default_settings_file(saved_state: Mapping[str, Any]) -> Path:
        save_file_location = saved_state.get("save_file_location")
        if isinstance(save_file_location, str) and save_file_location:
            return Path(save_file_location)
        return Path.cwd() / "settings.json"

    @staticmethod
    def write_file(file_path: str | Path, saved_state: Mapping[str, Any]) -> None:
        from persistence import write_project_settings_file, write_project_snapshot

        path = Path(file_path)
        if "ProjectArea" in saved_state and "ProjectSettings" in saved_state:
            project_settings = saved_state["ProjectSettings"]
            project_area = saved_state["ProjectArea"]
            if not isinstance(project_settings, object) or not isinstance(project_area, object):
                raise TypeError("Saved project state is missing native project objects.")
            session_state = saved_state.get("SessionState", {})
            form_state = saved_state.get("ProjectSettingsForm")
            write_project_snapshot(
                path,
                app_version=Activity.version,
                project_settings=project_settings,
                project_area=project_area,
                session_state=session_state if isinstance(session_state, Mapping) else {},
                form_state=form_state if isinstance(form_state, Mapping) else None,
            )
            return

        form_state = saved_state.get("ProjectSettingsForm")
        if isinstance(form_state, Mapping):
            write_project_settings_file(path, form_state, app_version=Activity.version)
            return

        raise ValueError("No saveable project settings or project snapshot data were available.")

    def no_clicked(self, checked: bool = False) -> None:
        del checked
        self.stop()


def start_save_state_activity(last_activity_type: type[Activity], saved_state: SavedState) -> None:
    del last_activity_type
    Activity.start_activity(SaveStateActivity, None, saved_state)
