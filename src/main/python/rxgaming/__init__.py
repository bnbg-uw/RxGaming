"""Tool for analyzing forest structure measured by lidar, comparing to reference conditions,
and gaming treatment scenarios.

The RxGaming tool is a graphical user interface for loading processed lidar data, summarizing
structure and pattern on the basis of individual trees, and comparing to a reference condition
dataset. Reference condition data can be accessed here (need to give URL).

The RxGaming tool works on the basis of stand polygons, optionally within a cohesive project area.
Structure is summarized for each stand and graphed against reference stands to assess departure.
Then, the possible range of treatments for each stand is graphed and the user can play out different
scenarios of treating for different residual structure and pattern."""

from . import activity
from . import QtWidgets
from . import projectsettings
from . import projectsettingsactivity
from . import gamingactivity

__all__ = ["activity", "QtWidgets", "projectsettings", "projectsettingsactivity", "gamingactivity"]
