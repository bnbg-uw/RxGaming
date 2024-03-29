from fbs_runtime.application_context.PyQt5 import ApplicationContext, cached_property
from activity import Activity, LoadStateActivity
from projectsettingsactivity import ProjectSettingsActivity
import sys
import traceback

class AppContext(ApplicationContext):           # 1. Subclass ApplicationContext
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


def handle_exception(exc_type, exc_value, exc_traceback):
    print("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
    sys.exit(1)


if __name__ == '__main__':
    sys.excepthook = handle_exception
    appctxt = AppContext()                      # 4. Instantiate the subclass
    prop_table = appctxt.get_resource("mcs_prop.csv")
    dll = appctxt.get_resource("bin/cpprxgaming.dll")
    appctxt.run(prop_table_path=prop_table, dll_path=dll)                   # 5. Invoke run()
    sys.exit()
