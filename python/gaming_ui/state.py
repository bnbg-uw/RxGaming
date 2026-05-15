from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping, Protocol, TypedDict, cast

from .units import UnitSystem


UnitSystemValue = Literal["metric", "imperial"]


class StandViewStatePayload(TypedDict):
    selected_unit_index: int
    active_page: int
    raster_mode: int
    show_treatment: bool
    dbh_cutoff: float
    unit_system: UnitSystemValue


def default_stand_view_state_payload() -> StandViewStatePayload:
    return {
        "selected_unit_index": 0,
        "active_page": 0,
        "raster_mode": 0,
        "show_treatment": False,
        "dbh_cutoff": 76.2,
        "unit_system": "imperial",
    }


def parse_stand_view_state_payload(raw: object, *, strict: bool = False) -> StandViewStatePayload:
    defaults = default_stand_view_state_payload()
    if raw is None:
        return defaults
    if not isinstance(raw, Mapping):
        if strict:
            raise ValueError("Session state must be a JSON object.")
        return defaults
    raw_mapping = cast(Mapping[str, object], raw)

    return {
        "selected_unit_index": _parse_int(raw_mapping, "selected_unit_index", defaults["selected_unit_index"], strict=strict),
        "active_page": _parse_int(raw_mapping, "active_page", defaults["active_page"], strict=strict),
        "raster_mode": _parse_int(raw_mapping, "raster_mode", defaults["raster_mode"], strict=strict),
        "show_treatment": _parse_bool(raw_mapping, "show_treatment", defaults["show_treatment"], strict=strict),
        "dbh_cutoff": _parse_float(raw_mapping, "dbh_cutoff", defaults["dbh_cutoff"], strict=strict),
        "unit_system": _parse_unit_system(raw_mapping, "unit_system", defaults["unit_system"], strict=strict),
    }


def _parse_int(raw: Mapping[str, object], key: str, default: int, *, strict: bool) -> int:
    value = raw.get(key, default)
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if strict:
        raise ValueError(f"Session state field '{key}' must be an integer.")
    return default


def _parse_bool(raw: Mapping[str, object], key: str, default: bool, *, strict: bool) -> bool:
    value = raw.get(key, default)
    if isinstance(value, bool):
        return value
    if strict:
        raise ValueError(f"Session state field '{key}' must be a boolean.")
    return default


def _parse_float(raw: Mapping[str, object], key: str, default: float, *, strict: bool) -> float:
    value = raw.get(key, default)
    if isinstance(value, int | float) and not isinstance(value, bool):
        return float(value)
    if strict:
        raise ValueError(f"Session state field '{key}' must be a number.")
    return default


def _parse_unit_system(
    raw: Mapping[str, object],
    key: str,
    default: UnitSystemValue,
    *,
    strict: bool,
) -> UnitSystemValue:
    value = raw.get(key, default)
    if value in ("metric", "imperial"):
        return value
    if strict:
        raise ValueError(f"Session state field '{key}' must be 'metric' or 'imperial'.")
    return default


@dataclass
class StandViewState:
    selected_unit_index: int = 0
    active_page: int = 0
    raster_mode: int = 0
    show_treatment: bool = False
    dbh_cutoff: float = 76.2
    unit_system: UnitSystem = UnitSystem.IMPERIAL

    def to_dict(self) -> StandViewStatePayload:
        return {
            "selected_unit_index": int(self.selected_unit_index),
            "active_page": int(self.active_page),
            "raster_mode": int(self.raster_mode),
            "show_treatment": bool(self.show_treatment),
            "dbh_cutoff": float(self.dbh_cutoff),
            "unit_system": "metric" if self.unit_system == UnitSystem.METRIC else "imperial",
        }

    @classmethod
    def from_dict(cls, payload: StandViewStatePayload) -> "StandViewState":
        return cls(
            selected_unit_index=payload["selected_unit_index"],
            active_page=payload["active_page"],
            raster_mode=payload["raster_mode"],
            show_treatment=payload["show_treatment"],
            dbh_cutoff=payload["dbh_cutoff"],
            unit_system=UnitSystem(payload["unit_system"]),
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
        return StandViewState.from_dict(parse_stand_view_state_payload(saved_state.get("SessionState")))

    def save_session(self, state: StandViewState, reason: str = "session_updated") -> None:
        del state
        del reason
