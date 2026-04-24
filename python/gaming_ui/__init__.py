from __future__ import annotations

from typing import Any

__all__ = ["GamingTabs"]


def __getattr__(name: str) -> Any:
    if name == "GamingTabs":
        from .tabs import GamingTabs

        return GamingTabs
    raise AttributeError(name)
