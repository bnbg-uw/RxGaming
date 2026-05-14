
"""
Copyright (C) 2024  University of Washington
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller & Sean Jeronimo
University of Washington Forest Resilience Lab
12/6/2024

activity.py
Custom framework for managing Qt windowing and application lifecycle, as well as saving/loading
persistent state information.
"""
from __future__ import annotations
from abc import ABC, abstractmethod
from enum import Enum, auto
from pathlib import Path
import sys
from typing import Any, Callable, ClassVar, Mapping

from PySide6.QtCore import QFileInfo, Qt
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication
from PySide6.QtWidgets import QFileIconProvider
from PySide6.QtWidgets import QMessageBox

from widgets import QWindow


SavedState = dict[str, Any]

def _resolve_application_icon() -> QIcon:
    icon_candidates = [
        Path(__file__).resolve().parents[1] / "icons" / "Icon.ico",
    ]

    if getattr(sys, "frozen", False):
        exe_dir = Path(sys.executable).resolve().parent
        icon_candidates.insert(0, exe_dir / "icons" / "Icon.ico")

    for candidate in icon_candidates:
        if candidate.exists():
            icon = QIcon(str(candidate))
            if not icon.isNull():
                return icon

    if sys.platform.startswith("win") and getattr(sys, "frozen", False):
        icon = QFileIconProvider().icon(QFileInfo(str(Path(sys.executable).resolve())))
        if not icon.isNull():
            return icon

    return QIcon()


def _get_qapplication() -> QApplication:
    instance = QApplication.instance()
    if instance is None:
        app = QApplication([])
    else:
        assert isinstance(instance, QApplication)
        app = instance
    icon = _resolve_application_icon()
    if not icon.isNull():
        app.setWindowIcon(icon)
    # The activity framework controls shutdown explicitly, so root-window
    # transitions should not make Qt quit underneath us.
    app.setQuitOnLastWindowClosed(False)
    return app


class WindowMode(Enum):
    """Enumeration of relationships between parent activity and child activities."""
    
    Modal = auto()
    """Parent window persists when child is created but is unresponsive until child window is killed."""
    Sibling = auto()
    """Child window's parent is set to parent window's parent. Both persist independently."""
    SimultaneousParent = auto()
    """Child and parent coexist but child is killed if parent is killed."""
    ExclusiveParent = auto() 
    """Parent is killed until child is killed, then parent is brought back.
    
    .. note:: With ``WindowMode.ExclusiveParent``, when the child is killed the parent that comes back is a
       newly created ``Activity`` re-instantiated from the ``saved_state`` of the old parent object.
    """

# Activity class (abstract base class)
#  Activities represent the basic unit of UI organization
#  Equivalent to one task, screen, view, etc.
#  Activities are created using the static start_activity factory method
#  Concrete subclasses need to implement on_start and save methods
class Activity(ABC):
    """Abstract base class for creating and managing activities.
    
    An activity is a single window with a simple lifecycle and the ability to save to/load from persistent states.
    Activities are defined by inheriting from this class and implementing the ``on_start`` and ``save`` methods.
    They are created via the ``Activity.start_activity`` static factory method, not by direct instantiation."""
    
    _app: ClassVar[QApplication] = _get_qapplication()
    _activities: ClassVar[list[Activity]] = []  # List of current activities. Exists so that activities don't get garbage collected
    _saved_state: ClassVar[SavedState] = {}  # Key/value pairs for all info that should be persistent across instances
    _starting: ClassVar[bool] = False
    _stopping: ClassVar[bool] = False  # Flag, True when application is shutting down
    _event_loop_running: ClassVar[bool] = False
    _last_activity_handler: ClassVar[Callable[[type[Activity], SavedState], None] | None] = None
    version = "1.0.11"
    
    try_to_save = False
    """Set to ``True`` if you want the application to prompt the user to save when exiting."""
    
    def __init__(self, parent_activity: Activity | None, window_mode: WindowMode):
        self._parent = parent_activity
        self._window_mode = window_mode
        self.window = QWindow()
        if not Activity._app.windowIcon().isNull():
            self.window.setWindowIcon(Activity._app.windowIcon())
        if self._window_mode is WindowMode.Modal:
            self.window.setWindowModality(Qt.WindowModality.ApplicationModal)  # TODO: set window hierarchy right so we can use WindowModal
        self.window.set_on_closed(self._on_window_close)

    def set_window(self, window: QWindow) -> None:
        self.window = window
        if not Activity._app.windowIcon().isNull():
            self.window.setWindowIcon(Activity._app.windowIcon())
        if self.window.on_closed is None:
            self.window.set_on_closed(self._on_window_close)
    
    # Each concrete subclass must implement this
    # Called after the window is set up
    # Responsible for setting up UI, callbacks, and loading from saved_state
    @abstractmethod
    def on_start(self, saved_state: SavedState, **kwargs: Any) -> None:  # TODO: think about passing window as argument
        """Lifecycle event callback, called after the `Activity`'s window is set up.
        
        Each concrete subclass must implement this method. This is the place to set up the UI,
        connect signals to slots for UI elements, and load from ``saved_state``.
        By the time this method is called, ``self.window`` is a fully functional, but not visible,
        ``QtWidgets.QWindow``. Setting up the UI usually goes something like this::
        
            def on_start(self, saved_state):
                self.label = QLabel("Don't push the button")
                self.button = QPushButton("Push me!")
                self.layout = QVBoxLayout()
                self.layout.addWidget(self.label)
                self.layout.addWidget(self.button)
                self.window.setLayout(self.layout)
                ...
        
        .. note:: UI elements must be class members so they are not garbage-collected after they go
            out of scope at the end of ``on_start``. Alternatively, they can be heap allocated, like this::
            
                self.layout.addWidget(QLabel("Don't push the button"))
                self.layout.addWidget(QPushButton("Push me!"))
        
        Next order of business is connecting signals to slots so the UI actually does something::
        
            def on_start(self, saved_state):
                ...
                self.button.clicked.connect(self.clicked)
                ...
            
            def clicked(self, event):
                self.label.setText("You shouldn't have!")
        
        And, lastly, loading from the ``saved_state`` dictionary::
        
            def on_start(self, saved_state):
                ...
                if "label_text" in saved_state:
                    label.setText(saved_state["label_text"])
        """
        ...
    
    @abstractmethod
    def save(self) -> SavedState:
        """Lifecycle event callback, called after the window is closed.
        
        Each concrete subclass must implement this. This is the time to save off all data required to
        recreate the current state. Ideally an activity that is stopped and saved off will be as similar
        as possible when re-instantiated.
        
        The saved state is to be returned as a ``dict``, where keys need to be known only by the ``on_start``
        method and values are whatever objects are required. For example, using the object described above::
        
            def save(self):
                return {"label_text": self.label.text()}
        """
        ...
    
    def stop(self) -> None:
        """Stops activity, closes window, and returns control to parent.
        
        If this activity has no parent, then control goes to the most recently created rootless
        (parent is ``None``) running activity. If there are no running activities, then the application
        will quit. In this event, if ``try_to_save`` is ``True``, then a ``SaveStateActivity`` will be
        created to handle saving off instances.
        """
        self.window.close()
    
    def _on_window_close(self, window: QWindow) -> None:
        """Callback for window closed event.
        
        Hooks in here to get instance, then passes to static method.
        """
        self._finish_activity(self)

    @classmethod
    def set_last_activity_handler(
        cls,
        handler: Callable[[type[Activity], SavedState], None] | None,
    ) -> None:
        """Registers project-specific behavior for final activity shutdown."""
        cls._last_activity_handler = handler

    @staticmethod
    def _handle_last_activity_closed(activity: Activity) -> None:
        '''Conceptually, this is created to allow hooking in to a saving on exit system, which may be project specific, so is not included in the ABC'''
        if Activity.try_to_save and Activity._last_activity_handler is not None:
            Activity._stopping = True
            Activity._saved_state.update({"LastActivity": type(activity)})
            Activity._last_activity_handler(type(activity), dict(Activity._saved_state))
            return
        Activity._app.quit()
    
    # Clean up after an activity
    #  Deal with saving state
    #  Remove from list of activities (thus it will be garbage-collected)
    #  Return focus to parent window
    @staticmethod
    def _finish_activity(activity: Activity) -> None:
        """Clean up after an activity.
        
        . Deal with saving state
        . Remove from list of activities (thus it will be garbage-collected)
        . Return focus to parent window or quit application
        """
        Activity._saved_state.update(activity.save())
        if activity in Activity._activities:
            Activity._activities.remove(activity)
        if Activity._stopping:
            if len(Activity._activities) == 0:
                Activity._stopping = False
                Activity._app.quit()
            return
        if Activity._starting or Activity._stopping:
            return
        
        # Deal with child activities -- kill off children, reassign siblings
        children = [a for a in Activity._activities if a._parent is activity]
        for c in children:
            if c._window_mode is WindowMode.Sibling:
                c._parent = activity._parent
            if c._window_mode is WindowMode.SimultaneousParent or c._window_mode is WindowMode.ExclusiveParent:
                c.stop()
        
        # Figure out who gets control next
        if activity._parent is None:  # In this case, find the last activity in the list that is a root activity
            root_activities = [a for a in Activity._activities if a._parent is None]
            if len(root_activities) > 0:
                root_activities[-1].window.activateWindow()
        elif activity._window_mode is WindowMode.Modal or activity._window_mode is WindowMode.Sibling or \
             activity._window_mode is WindowMode.SimultaneousParent:
            activity._parent.window.activateWindow()
        elif activity._window_mode is WindowMode.ExclusiveParent:
            Activity.start_activity(type(activity._parent), activity._parent._parent, Activity._saved_state)
        
        # If all activities are closed, either save and quit or just stop the
        # current event loop so the caller can decide what to do next.
        if len(Activity._activities) == 0:
            Activity._handle_last_activity_closed(activity)
    
    # Factory method for starting a new activity
    # Call this instead of creating Activity instances directly through constructor
    # Pass in a parent activity unless the Activity you are starting will be the root
    @staticmethod
    def start_activity(
        activity_class_to_start: type[Activity],
        parent_activity: Activity | None = None,
        saved_state: Mapping[str, Any] | None = None,
        window_mode: WindowMode = WindowMode.SimultaneousParent,
        **kwargs: Any,
    ) -> Activity | None:
        """Factory method for creating and starting a new activity.
        
        Call this instead of creating ``Activity`` instances through the constructor.
        If you do not override the default for ``parent_activity``, then the created activity
        will be a root activity with no specific relationship to other activities/windows.
        You may pass in a ``saved_state``, which is usually for reconstructing states saved by
        ``SaveStateActivity`` in a previous session, but can also be used to pass in arguments
        to the new activity's ``on_start`` method.
        
        .. note:: This method starts and handles the Qt event loop. The event loop does not need to,
           and should not, be started again anywhere else.

        See documentation for |WindowMode link|_ for information on the possible relationships between
        parent and child activities.
        
        .. |WindowMode link| replace:: ``WindowMode``
        .. _WindowMode link: activity.WindowMode.html
        """
        Activity._starting = True
        state_update = dict(saved_state or {})
        Activity._saved_state.update(state_update)
        current_activity = activity_class_to_start(parent_activity, window_mode)

        try:
            current_activity.on_start(dict(Activity._saved_state), **kwargs)
            if window_mode is WindowMode.ExclusiveParent and parent_activity is not None:
                parent_activity.stop()
            current_activity.window.show()
            current_activity.window.activateWindow()
            Activity._activities.append(current_activity)
        finally:
            Activity._starting = False

        if not Activity._event_loop_running:
            if hasattr(Activity, "test"):
                return current_activity
            Activity._event_loop_running = True
            try:
                Activity._app.exec()
            finally:
                Activity._event_loop_running = False

        return current_activity

    @staticmethod
    def _show_message(
        icon: QMessageBox.Icon,
        title: str,
        text: str,
        informative_text: str | None = None,
    ) -> None:
        msg = QMessageBox()
        msg.setIcon(icon)
        msg.setText(text)
        if informative_text is not None:
            msg.setInformativeText(informative_text)
        msg.setWindowTitle(title)
        if not Activity._app.windowIcon().isNull():
            msg.setWindowIcon(Activity._app.windowIcon())
        msg.setStandardButtons(QMessageBox.StandardButton.Ok)
        msg.exec()
