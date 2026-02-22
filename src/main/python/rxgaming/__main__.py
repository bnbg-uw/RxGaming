"""
    Copyright (C) 2024  University of Washington
    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see https://www.gnu.org/licenses/.

Bryce Bartl-Geller
University of Washington Forest Resilience Lab
12/6/2024

__main__.py
This is the entry point to the tool. We use fman build system to build the executable: https://build-system.fman.io/
"""

from fbs_runtime.application_context.PyQt5 import ApplicationContext, cached_property
from activity import Activity, LoadStateActivity
from projectsettingsactivity import ProjectSettingsActivity
import sys
import traceback

# This makes fman build system work with our Acitivty class which is how we switch between windows in this program.
class AppContext(ApplicationContext):
    def __init__(self):
        super(AppContext, self).__init__()
        self._to_start = ProjectSettingsActivity

    def run(self, **kwargs):                              # 2. Implement run()
        Activity.Start_Activity(LoadStateActivity, saved_state={'onLoad': self.onLoad})
        Activity.Try_To_Save = True
        Activity.Start_Activity(self._to_start, **kwargs)


    def onLoad(self, saved_state):
        if "LastActivity" in saved_state:
            self._to_start = saved_state["LastActivity"]

    @cached_property
    def app(self):
        return Activity._App


# In theory we need something more elegant here.
def handle_exception(exc_type, exc_value, exc_traceback):
    print("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    sys.exit(1)


if __name__ == '__main__':
    sys.excepthook = handle_exception
    appctxt = AppContext()
    prop_table = appctxt.get_resource("mcs_prop.csv")
    dll = appctxt.get_resource("bin/pyRxTools.dll")
    appctxt.run(prop_table_path=prop_table, dll_path=dll)
    sys.exit()
