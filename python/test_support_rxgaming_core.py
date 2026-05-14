from __future__ import annotations

from pathlib import Path
import sys
import types


class StubProjectSettings:
    def __init__(self, *args: object, **kwargs: object) -> None:
        fields = [
            "name",
            "unitPolyPath",
            "refDataPath",
            "mcsPropPath",
            "fiaPath",
            "lidarPath",
            "unitName",
            "savePath",
            "nThread",
        ]
        for field, value in zip(fields, args):
            setattr(self, field, value)
        for field in fields:
            if field in kwargs:
                setattr(self, field, kwargs[field])


class StubProjectArea:
    def __init__(self, *args: object, **kwargs: object) -> None:
        del args, kwargs


class StubProjectAreaBuildSnapshot:
    def __init__(self) -> None:
        self.stage = ""
        self.message = ""
        self.completed = -1
        self.total = -1
        self.status = "running"
        self.error = ""


def ensure_rxgaming_core_test_module() -> types.ModuleType:
    module = sys.modules.get("rxgaming_core")
    if module is None:
        try:
            import rxgaming_core as module  # type: ignore
        except ImportError:
            module = types.ModuleType("rxgaming_core")
            sys.modules["rxgaming_core"] = module

    module.ProjectSettings = getattr(module, "ProjectSettings", StubProjectSettings)
    module.ProjectArea = getattr(module, "ProjectArea", StubProjectArea)
    module.RxUnit = getattr(module, "RxUnit", object)
    module.StructureSummary = getattr(module, "StructureSummary", object)
    module.TreatmentEngine = getattr(module, "TreatmentEngine", object)
    module.ProjectAreaBuildHandle = getattr(module, "ProjectAreaBuildHandle", object)
    module.ProjectAreaBuildSnapshot = getattr(module, "ProjectAreaBuildSnapshot", StubProjectAreaBuildSnapshot)
    module.start_project_area_build = getattr(module, "start_project_area_build", lambda project_settings: object())
    module.poll_project_area_build = getattr(module, "poll_project_area_build", lambda handle: StubProjectAreaBuildSnapshot())
    module.finish_project_area_build = getattr(module, "finish_project_area_build", lambda handle: StubProjectArea())
    module.load_project_area = getattr(module, "load_project_area", lambda path: StubProjectArea())
    module.save_project_area = getattr(module, "save_project_area", lambda project_area, path: Path(path).write_bytes(b"stub"))
    module.set_proj_db_path = getattr(module, "set_proj_db_path", lambda path: None)
    return module
