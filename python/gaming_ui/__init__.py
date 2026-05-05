from __future__ import annotations

__all__ = ["GamingTabs"]


def __getattr__(name: str):
    if name == "GamingTabs":
        from .tabs import GamingTabs

        return GamingTabs
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
