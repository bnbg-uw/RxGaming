from __future__ import annotations

from enum import StrEnum
from typing import Any, Literal

STRUCTURE_METRICS: tuple[MetricKind, ...] = ("tph", "ba", "mcs", "cc")

# These constants intentionally mirror the values in src/lapisgis/src/Unit.hpp.
METER_PER_FOOT = 0.3048
METER_PER_INCH = 0.0254
METER_PER_CENTIMETER = 0.01
SQUARE_METER_PER_HECTARE = 10000.0
SQUARE_METER_PER_ACRE = 4046.8564224
PERCENT_PER_PROPORTION = 100.0

CENTIMETERS_PER_INCH = METER_PER_INCH / METER_PER_CENTIMETER
FEET_PER_METER = 1.0 / METER_PER_FOOT
TPA_PER_TPH = SQUARE_METER_PER_ACRE / SQUARE_METER_PER_HECTARE
FT2_PER_ACRE_PER_M2_PER_HA = FEET_PER_METER * FEET_PER_METER * TPA_PER_TPH

MetricKind = Literal["dbh", "height", "tph", "ba", "mcs", "cc", "canopy_height"]


class UnitSystem(StrEnum):
    METRIC = "metric"
    IMPERIAL = "imperial"


def area_to_display(metric_area_ha: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return metric_area_ha / TPA_PER_TPH
    return metric_area_ha


def area_label(unit_system: UnitSystem) -> str:
    return "ac" if unit_system == UnitSystem.IMPERIAL else "ha"


def format_area(metric_area_ha: float, unit_system: UnitSystem, precision: int = 2) -> str:
    display_value = area_to_display(metric_area_ha, unit_system)
    return f"{display_value:.{precision}f} {area_label(unit_system)}"


def dbh_to_display(metric_cm: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return metric_cm / CENTIMETERS_PER_INCH
    return metric_cm


def dbh_from_display(display_value: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return display_value * CENTIMETERS_PER_INCH
    return display_value


def height_to_display(metric_m: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return metric_m * FEET_PER_METER
    return metric_m


def height_from_display(display_value: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return display_value / FEET_PER_METER
    return display_value


def tph_to_display(metric_tph: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return metric_tph * TPA_PER_TPH
    return metric_tph


def tph_from_display(display_value: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return display_value / TPA_PER_TPH
    return display_value


def ba_to_display(metric_m2_per_ha: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return metric_m2_per_ha * FT2_PER_ACRE_PER_M2_PER_HA
    return metric_m2_per_ha


def ba_from_display(display_value: float, unit_system: UnitSystem) -> float:
    if unit_system == UnitSystem.IMPERIAL:
        return display_value / FT2_PER_ACRE_PER_M2_PER_HA
    return display_value


def label_for(metric_kind: MetricKind, unit_system: UnitSystem) -> str:
    if metric_kind == "dbh":
        return "in" if unit_system == UnitSystem.IMPERIAL else "cm"
    if metric_kind in {"height", "canopy_height"}:
        return "ft" if unit_system == UnitSystem.IMPERIAL else "m"
    if metric_kind == "tph":
        return "TPA" if unit_system == UnitSystem.IMPERIAL else "TPH"
    if metric_kind == "ba":
        return "ft^2/ac" if unit_system == UnitSystem.IMPERIAL else "m^2/ha"
    if metric_kind == "mcs":
        return "MCS"
    if metric_kind == "cc":
        return "%"
    raise ValueError(f"Unknown metric kind: {metric_kind}")


def display_name_for(metric_kind: MetricKind, unit_system: UnitSystem) -> str:
    if metric_kind == "dbh":
        return f"DBH ({label_for(metric_kind, unit_system)})"
    if metric_kind == "height":
        return f"Height ({label_for(metric_kind, unit_system)})"
    if metric_kind == "canopy_height":
        return f"Canopy Height Model ({label_for(metric_kind, unit_system)})"
    if metric_kind == "tph":
        return f"Density ({label_for(metric_kind, unit_system)})"
    if metric_kind == "ba":
        return f"Basal Area ({label_for(metric_kind, unit_system)})"
    if metric_kind == "mcs":
        return "Mean Clump Size"
    if metric_kind == "cc":
        return "Canopy Cover"
    raise ValueError(f"Unknown metric kind: {metric_kind}")


def to_display(metric_kind: MetricKind, metric_value: float, unit_system: UnitSystem) -> float:
    if metric_kind == "dbh":
        return dbh_to_display(metric_value, unit_system)
    if metric_kind in {"height", "canopy_height"}:
        return height_to_display(metric_value, unit_system)
    if metric_kind == "tph":
        return tph_to_display(metric_value, unit_system)
    if metric_kind == "ba":
        return ba_to_display(metric_value, unit_system)
    if metric_kind == "mcs":
        return metric_value
    if metric_kind == "cc":
        return metric_value * PERCENT_PER_PROPORTION
    raise ValueError(f"Unknown metric kind: {metric_kind}")


def from_display(metric_kind: MetricKind, display_value: float, unit_system: UnitSystem) -> float:
    if metric_kind == "dbh":
        return dbh_from_display(display_value, unit_system)
    if metric_kind in {"height", "canopy_height"}:
        return height_from_display(display_value, unit_system)
    if metric_kind == "tph":
        return tph_from_display(display_value, unit_system)
    if metric_kind == "ba":
        return ba_from_display(display_value, unit_system)
    if metric_kind == "mcs":
        return display_value
    if metric_kind == "cc":
        return display_value / PERCENT_PER_PROPORTION
    raise ValueError(f"Unknown metric kind: {metric_kind}")


def format_value(
    metric_kind: MetricKind,
    metric_value: float,
    unit_system: UnitSystem,
    precision: int = 2,
    include_unit: bool = False,
) -> str:
    display_value = to_display(metric_kind, metric_value, unit_system)
    formatted = f"{display_value:.{precision}f}"
    if include_unit:
        return f"{formatted} {label_for(metric_kind, unit_system)}"
    return formatted


def array_to_display(metric_kind: MetricKind, values: Any, unit_system: UnitSystem) -> Any:
    if metric_kind == "dbh":
        factor = 1.0 / CENTIMETERS_PER_INCH if unit_system == UnitSystem.IMPERIAL else 1.0
    elif metric_kind in {"height", "canopy_height"}:
        factor = FEET_PER_METER if unit_system == UnitSystem.IMPERIAL else 1.0
    elif metric_kind == "tph":
        factor = TPA_PER_TPH if unit_system == UnitSystem.IMPERIAL else 1.0
    elif metric_kind == "ba":
        factor = FT2_PER_ACRE_PER_M2_PER_HA if unit_system == UnitSystem.IMPERIAL else 1.0
    elif metric_kind == "mcs":
        factor = 1.0
    elif metric_kind == "cc":
        factor = PERCENT_PER_PROPORTION
    else:
        raise ValueError(f"Unknown metric kind: {metric_kind}")
    return values * factor
