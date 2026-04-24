from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Protocol

from .units import UnitSystem


@dataclass
class StandViewState:
    selected_unit_index: int = 0
    active_page: int = 0
    raster_mode: int = 0
    show_treatment: bool = False
    dbh_cutoff: float = 76.2
    unit_system: UnitSystem = UnitSystem.IMPERIAL

    def to_dict(self) -> dict[str, object]:
        return {
            "selected_unit_index": int(self.selected_unit_index),
            "active_page": int(self.active_page),
            "raster_mode": int(self.raster_mode),
            "show_treatment": bool(self.show_treatment),
            "dbh_cutoff": float(self.dbh_cutoff),
            "unit_system": self.unit_system.value,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> "StandViewState":
        defaults = cls()
        return cls(
            selected_unit_index=int(payload.get("selected_unit_index", 0)),
            active_page=int(payload.get("active_page", 0)),
            raster_mode=int(payload.get("raster_mode", 0)),
            show_treatment=bool(payload.get("show_treatment", False)),
            dbh_cutoff=float(payload.get("dbh_cutoff", defaults.dbh_cutoff)),
            unit_system=UnitSystem(payload.get("unit_system", defaults.unit_system.value)),
        )

    @classmethod
    def from_saved_state(cls, saved_state: dict[str, object]) -> "StandViewState":
        session_state = saved_state.get("SessionState")
        if isinstance(session_state, Mapping):
            return cls.from_dict(session_state)

        defaults = cls()
        cutoff = saved_state.get(
            "GamingActivity.cut_range",
            saved_state.get("GamingActivity.dbh_max", defaults.dbh_cutoff),
        )
        return cls(
            selected_unit_index=int(saved_state.get("GamingActivity.selected_unit", 0)),
            active_page=int(saved_state.get("GamingActivity.page", 0)),
            raster_mode=int(saved_state.get("GamingActivity.raster_mode", 0)),
            show_treatment=bool(saved_state.get("GamingActivity.show_treatment", False)),
            dbh_cutoff=float(cutoff),
            unit_system=UnitSystem(saved_state.get("GamingActivity.unit_system", defaults.unit_system.value)),
        )

    @property
    def dbh_min(self) -> float:
        return 0.0

    @property
    def dbh_max(self) -> float:
        return self.dbh_cutoff


class GamingSessionPersistence(Protocol):
    def load_initial_state(self, saved_state: dict[str, object]) -> StandViewState: ...

    def save_session(self, state: StandViewState, reason: str = "session_updated") -> None: ...


class NoOpGamingSessionPersistence:
    def load_initial_state(self, saved_state: dict[str, object]) -> StandViewState:
        return StandViewState.from_saved_state(saved_state)

    def save_session(self, state: StandViewState, reason: str = "session_updated") -> None:
        del state
        del reason
