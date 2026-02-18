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

from enum import Enum, auto
from abc import *
import pickle
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout, QLabel, QPushButton, QFileDialog, QMessageBox
from PyQt5.QtCore import Qt
from QtWidgets import QWindow

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


# Define dummy class to use for function signature type specs
class Activity:
    pass

# Activity class (abstract base class)
#  Activities represent the basic unit of UI organization
#  Equivalent to one task, screen, view, etc.
#  Activities are created using the static Start_Activity factory method
#  Concrete subclasses need to implement on_start and save methods
class Activity(ABC):
    """Abstract base class for creating and managing activities.
    
    An activity is a single window with a simple lifecycle and the ability to save to/load from persistent states.
    Activities are defined by inheriting from this class and implementing the ``on_start`` and ``save`` methods.
    They are created via the ``Activity.Start_Activity`` static factory method, not by direct instantiation."""
    
    _App = QApplication([])  # The single instance of Qt application
    _Activities = []  # List of current activities. Exists so that activities don't get garbage collected
    _Saved_State = {}  # Key/value pairs for all info that should be persistent across instances
    _Stopping = False  # Flag, True when application is shutting down
    version = "1.0.11"
    
    Try_To_Save = False
    """Set to ``True`` if you want the application to prompt the user to save when exiting."""
    
    def __init__(self, parent_activity: Activity, window_mode: WindowMode):
        self._parent = parent_activity
        self._window_mode = window_mode
        self.window = QWindow()
        if self._window_mode is WindowMode.Modal:
            self.window.setWindowModality(Qt.ApplicationModal)  # TODO: set window hierarchy right so we can use WindowModal
        self.window.setOnClosed(self._onWindowClose)

    def set_window(self, window):
        self.window = window
        if self.window.onClosed is None:
            self.window.setOnClosed(self._onWindowClose)
    
    # Each concrete subclass must implement this
    # Called after the window is set up
    # Responsible for setting up UI, callbacks, and loading from saved_state
    @abstractmethod
    def on_start(self, saved_state: dict):  # TODO: think about passing window as argument
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
    def save(self) -> dict:
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
    
    def stop(self):
        """Stops activity, closes window, and returns control to parent.
        
        If this activity has no parent, then control goes to the most recently created rootless
        (parent is ``None``) running activity. If there are no running activities, then the application
        will quit. In this event, if ``Try_To_Save`` is ``True``, then a ``SaveStateActivity`` will be
        created to handle saving off instances.
        """
        self.window.close()
    
    def _onWindowClose(self, window: QWindow):
        """Callback for window closed event.
        
        Hooks in here to get instance, then passes to static method.
        """
        Activity._FinishActivity(self)
    
    # Clean up after an activity
    #  Deal with saving state
    #  Remove from list of activities (thus it will be garbage-collected)
    #  Return focus to parent window
    @staticmethod
    def _FinishActivity(activity: Activity):
        """Clean up after an activity.
        
        . Deal with saving state
        . Remove from list of activities (thus it will be garbage-collected)
        . Return focus to parent window or quit application
        """
        Activity._Saved_State.update(activity.save())
        Activity._Activities.remove(activity)
        if Activity._Starting or Activity._Stopping:
            return
        
        # Deal with child activities -- kill off children, reassign siblings
        children = [a for a in Activity._Activities if a._parent is activity]
        for c in children:
            if c._window_mode is WindowMode.Sibling:
                c._parent = activity._parent
            if c._window_mode is WindowMode.SimultaneousParent or c._window_mode is WindowMode.ExclusiveParent:
                c.stop()
        
        # Figure out who gets control next
        if activity._parent is None:  # In this case, find the last activity in the list that is a root activity
            root_activities = [a for a in Activity._Activities if a._parent is None]
            if len(root_activities) > 0:
                root_activities[-1].window.activateWindow()
        elif activity._window_mode is WindowMode.Modal or activity._window_mode is WindowMode.Sibling or \
             activity._window_mode is WindowMode.SimultaneousParent:
            activity._parent.window.activateWindow()
        elif activity._window_mode is WindowMode.ExclusiveParent:
            Activity.Start_Activity(type(activity._parent), activity._parent._parent, Activity._Saved_State)
        
        # If all activities are closed, then we save and quit
        if len(Activity._Activities) == 0 and Activity.Try_To_Save:
            Activity._Stopping = True
            Activity._Saved_State.update({"LastActivity": type(activity)})
            Activity.Start_Activity(SaveStateActivity, None, Activity._Saved_State)
    
    # Factory method for starting a new activity
    # Call this instead of creating Activity instances directly through constructor
    # Pass in a parent activity unless the Activity you are starting will be the root
    @staticmethod
    def Start_Activity(activity_class_to_start: Activity, parent_activity: Activity=None,
                       saved_state: dict={}, window_mode: WindowMode=WindowMode.SimultaneousParent, **kwargs):
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
        Activity._Starting = True
        Activity._Saved_State.update(saved_state)
        current_activity = activity_class_to_start(parent_activity, window_mode)
        current_activity.on_start(Activity._Saved_State, **kwargs)
        if window_mode is WindowMode.ExclusiveParent and parent_activity is not None:
            parent_activity.stop()
        current_activity.window.show()
        current_activity.window.activateWindow()
        Activity._Activities.append(current_activity)
        Activity._Starting = False
        if Activity._App.applicationState() == Qt.ApplicationInactive:
            if hasattr(Activity, "test"):
                return current_activity
            else:
                Activity._App.exec_()


class SaveStateActivity(Activity):
    """Prompts user to store the state instance dictionary.
    
    The ``SaveStateActivity`` generally does not need to be created,
    since it is primarily used internally by ``Activity``.
    However, there is no reason it couldn't be used if desired.
    
    Usage is like any other activity::
    
        Activity.Start_Activity(SaveStateActivity)
    """
    
    def on_start(self, saved_state):
        """"""
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
    
    def save(self):
        """"""
        return {}
    
    def yes_clicked(self, button):
        self.prompt_and_save(self.saved_state)
        self.stop()

    def prompt_and_save(self, saved_state):
        file_path = QFileDialog.getSaveFileName(None, "Save as...", "", "*.dat")[0]
        if file_path == "":
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Please select a file (.dat) to save your work")
            msg.setWindowTitle("Choose a file")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            return
        self.write_file(file_path, saved_state)

    @staticmethod
    def write_file(file_path, saved_state):
        with open(file_path, "wb") as fp:
            pickle.dump(saved_state, fp)

    def no_clicked(self, button):
        self.stop()


# LoadStateActivity class
#  Prompts user to begin new project or load from saved state
#  If user chooses to load, it populates the Activity state instance
class LoadStateActivity(Activity):
    """Prompts user to begin new project or load from saved state.
    
    If user chooses to start a new project, this activity quits and it is assumed that the desired
    activity will be started next. If the user chooses to load, then a file dialog is shown,
    the selected file is loaded, and the loaded state instance is made available for subsequent activities.
    
    This activity looks for a special item in the ``saved_state`` dictionary (`onLoad`) that can provide an
    optional callback to be executed once the state has been loaded. This can be useful when deciding what
    part of an application to start given the state upon exit of the last session. There is also a special
    item, `LastActivity`, that is saved whenever an application exits and can be used for this end.
    Here is an example using this model::
    
        def on_loaded(saved_state):
            if "LastActivity" in saved_state:
                main.to_start = saved_state["Last_Activity"]
        
        def main():
            main.to_start = NewProjectActivity
            Activity.Start_Activity(LoadStateActivity, saved_state={"onLoad": on_loaded})
            Activity.Start_Activity(main.to_start)
    """     
    
    def on_start(self, saved_state):
        """"""
        self.onLoad = saved_state["onLoad"] if "onLoad" in saved_state else None
        Activity._Saved_State = {}
        
        # TODO make strings a declarative setting/localizable
        self.label = QLabel("What would you like to do?")
        self.load_button = QPushButton("Work from saved file")
        self.new_button = QPushButton("Start new project")
        
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.load_button)
        self.layout.addWidget(self.new_button)
        self.window.setLayout(self.layout)
        self.window.setWindowTitle("Load a file")
        
        self.load_button.clicked.connect(self.load_clicked)
        self.new_button.clicked.connect(self.new_clicked)
    
    def save(self):
        """"""
        return {}
    
    def load_clicked(self, button):
        file_path = QFileDialog.getOpenFileName(self.window, "Open...", "", "*.dat")[0]
        if file_path == "":
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Warning)
            msg.setText("Please select a file (.dat) to open")
            msg.setWindowTitle("Choose a file")
            msg.setStandardButtons(QMessageBox.Ok)
            msg.exec_()
            return
        with open(file_path, "rb") as fp:
            saved_state = pickle.load(fp)
            saved_state['save_file_location'] = file_path
        Activity._Saved_State = saved_state
        if self.onLoad is not None:
            self.onLoad(saved_state)
        self.stop()
    
    def new_clicked(self, button):
        self.stop()
