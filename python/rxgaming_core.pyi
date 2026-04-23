from __future__ import annotations

from enum import Enum
from typing import Sequence

import numpy as np
import numpy.typing as npt


FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.integer]


def set_proj_db_path(path: str) -> None: ...
def set_seed(seed: int) -> None: ...
def get_ba_dist(taos: FloatArray) -> FloatArray: ...


class StructureSummary:
    ba: float
    tph: float
    mcs: float
    cc: float
    csd: Sequence[float]
    binMins: Sequence[float]
    binMaxs: Sequence[float]

    def __getitem__(self, i: int) -> float: ...


class TreatmentResult(Enum):
    success: TreatmentResult
    diameterFailure: TreatmentResult
    cuttingFailure: TreatmentResult


class ProjectSettings:
    name: str
    unitPolyPath: str
    refDataPath: str
    mcsPropPath: str
    fiaPath: str
    lidarPath: str
    unitName: str
    savePath: str
    nThread: int

    def __init__(
        self,
        name: str,
        unitPolyPath: str,
        refDataPath: str,
        mcsPropPath: str,
        fiaPath: str,
        lidarPath: str,
        unitName: str,
        savePath: str,
        nThread: int,
    ) -> None: ...


class RxUnit:
    name: str
    areaHa: float
    result: TreatmentResult
    currentStructure: StructureSummary
    targetStructure: StructureSummary
    treatedStructure: StructureSummary

    def get_mask(self) -> IntArray: ...
    def get_chm(self) -> FloatArray: ...
    def get_basin(self) -> IntArray: ...
    def get_hillshade(self) -> FloatArray: ...
    def get_clump_map(self) -> IntArray: ...
    def get_taos(self) -> FloatArray: ...
    def get_clump_sizes(self) -> FloatArray: ...
    def get_treat_chm(self) -> FloatArray: ...
    def get_treat_basin(self) -> IntArray: ...
    def get_treat_hillshade(self) -> FloatArray: ...
    def get_treat_clump_map(self) -> IntArray: ...
    def get_treat_taos(self) -> FloatArray: ...
    def get_treat_clump_sizes(self) -> FloatArray: ...
    def get_cut_taos(self) -> FloatArray: ...
    def get_simulated_structures(self, bbDbh: float) -> list[StructureSummary]: ...


class TreatmentEngine:
    def __init__(self) -> None: ...
    def do_treatment(self, unit: RxUnit, dbhMin: float, dbhMax: float) -> None: ...


class ProjectArea:
    rxUnits: list[RxUnit]

    def __init__(self, ps: ProjectSettings) -> None: ...
