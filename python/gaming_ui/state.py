from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from .units import UnitSystem


@dataclass(slots=True)
class StandViewState:
    selected_unit_index: int = 0
    active_page: int = 0
    raster_mode: int = 0
    show_treatment: bool = False
    dbh_cutoff: float = 76.2
    preview_dbh: float = 76.2
    unit_system: UnitSystem = UnitSystem.IMPERIAL

    @classmethod
    def from_saved_state(cls, saved_state: dict[str, object]) -> "StandViewState":
        cutoff = saved_state.get("GamingActivity.cut_range", saved_state.get("GamingActivity.dbh_max", 30.0))
        return cls(
            selected_unit_index=int(saved_state.get("GamingActivity.selected_unit", 0)),
            active_page=int(saved_state.get("GamingActivity.page", 0)),
            raster_mode=int(saved_state.get("GamingActivity.raster_mode", 0)),
            show_treatment=bool(saved_state.get("GamingActivity.show_treatment", False)),
            dbh_cutoff=float(cutoff),
            preview_dbh=float(saved_state.get("GamingActivity.preview_dbh", 30.0)),
            unit_system=UnitSystem(saved_state.get("GamingActivity.unit_system", UnitSystem.IMPERIAL.value)),
        )

    @property
    def dbh_min(self) -> float:
        return 0.0

    @property
    def dbh_max(self) -> float:
        return self.dbh_cutoff


class GamingSessionPersistence(Protocol):
    def load_initial_state(self, saved_state: dict[str, object]) -> StandViewState: ...

    def save_session(self, state: StandViewState) -> None: ...


class NoOpGamingSessionPersistence:
    def load_initial_state(self, saved_state: dict[str, object]) -> StandViewState:
        return StandViewState.from_saved_state(saved_state)

    def save_session(self, state: StandViewState) -> None:
        del state
