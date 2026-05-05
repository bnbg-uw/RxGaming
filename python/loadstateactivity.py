from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtWidgets import QFileDialog, QLabel, QMessageBox, QPushButton, QVBoxLayout

from activity import Activity, SavedState


class LoadStateActivity(Activity):
    """Prompts user to begin new project or load from saved state."""

    def on_start(self, saved_state: SavedState, **kwargs: Any) -> None:
        del kwargs
        self.on_load = saved_state.get("onLoad")
        self._continue_after_close = False
        Activity._saved_state = {}

        self.label = QLabel("What would you like to do?")
        self.load_project_button = QPushButton("Work from saved project")
        self.load_settings_button = QPushButton("Load project settings")
        self.new_button = QPushButton("Start new project")

        self.layout = QVBoxLayout()
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.load_project_button)
        self.layout.addWidget(self.load_settings_button)
        self.layout.addWidget(self.new_button)
        self.window.setLayout(self.layout)
        self.window.setWindowTitle("Load a file")

        self.load_project_button.clicked.connect(self.load_project_clicked)
        self.load_settings_button.clicked.connect(self.load_settings_clicked)
        self.new_button.clicked.connect(self.new_clicked)

    def save(self) -> SavedState:
        return {"LoadStateContinue": self._continue_after_close}

    def load_project_clicked(self, checked: bool = False) -> None:
        del checked
        selected_path = QFileDialog.getExistingDirectory(self.window, "Select project folder")
        if selected_path == "":
            Activity._show_message(
                QMessageBox.Icon.Warning,
                "Choose a location",
                "Please select a project folder to open",
            )
            return
        path = Path(selected_path)

        from persistence import read_project_snapshot

        try:
            loaded = read_project_snapshot(path)
            saved_state = {
                "ProjectSettings": loaded.project_settings,
                "ProjectArea": loaded.project_area,
                "ProjectSettingsForm": loaded.form_state,
                "ProjectSnapshotPath": str(loaded.project_root),
                "SessionState": loaded.session_state,
            }
        except Exception as exc:
            Activity._show_message(
                QMessageBox.Icon.Warning,
                "Open failed",
                str(exc),
            )
            return

        Activity._saved_state = saved_state
        self._continue_after_close = True
        if self.on_load is not None:
            self.on_load(saved_state)
        self.stop()

    def load_settings_clicked(self, checked: bool = False) -> None:
        del checked
        selected_path = QFileDialog.getOpenFileName(self.window, "Select settings .json", "", "*.json")[0]
        if selected_path == "":
            Activity._show_message(
                QMessageBox.Icon.Warning,
                "Choose a location",
                "Please select a settings JSON file to open to open",
            )
            return
        path = Path(selected_path)

        if path.is_file() and path.suffix.lower() != ".json":
            Activity._show_message(
                QMessageBox.Icon.Warning,
                "Unsupported file",
                "Please choose a settings JSON file or a project folder.",
            )
            return

        from persistence import read_project_settings_file

        try:
            loaded = read_project_settings_file(path)
            saved_state = {
                "ProjectSettingsForm": loaded.form_state,
                "save_file_location": str(loaded.settings_path),
            }
        except Exception as exc:
            Activity._show_message(
                QMessageBox.Icon.Warning,
                "Open failed",
                str(exc),
            )
            return

        Activity._saved_state = saved_state
        self._continue_after_close = True
        if self.on_load is not None:
            self.on_load(saved_state)
        self.stop()

    def new_clicked(self, checked: bool = False) -> None:
        del checked
        self._continue_after_close = True
        self.stop()
