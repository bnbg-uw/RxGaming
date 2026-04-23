from __future__ import annotations

from typing import Any

import numpy as np
from matplotlib.figure import Figure
from PySide6.QtWidgets import QTabWidget, QWidget

from rxgaming_core import ProjectArea, ProjectSettings, RxUnit
from .landscape import LandscapeReferenceTab
from .stand import StandViewCoordinator
from .state import GamingSessionPersistence, NoOpGamingSessionPersistence
from .units import UnitSystem


class GamingTabs(QTabWidget):
    def __init__(
        self,
        project_settings: ProjectSettings,
        project_area: ProjectArea,
        saved_state: dict[str, Any],
        persistence: GamingSessionPersistence | None = None,
    ) -> None:
        super().__init__(None)
        self.rx_units = project_area.rxUnits
        self._persistence = persistence or NoOpGamingSessionPersistence()
        self._state = self._persistence.load_initial_state(saved_state)

        self.stand_tab = StandViewCoordinator(self.rx_units, self._state)
        self.landscape_tab = LandscapeReferenceTab(self.rx_units, project_settings, self._state)

        self.stand_tab.set_state_changed_callback(self._handle_stand_state_change)
        self.landscape_tab.set_unit_selected_callback(self.stand_tab.select_unit)

        self.addTab(self.stand_tab, "Stand View")
        self.addTab(self.landscape_tab, "Landscape View")
        self._handle_stand_state_change()

    @property
    def raster_figure(self) -> Figure:
        return self.stand_tab.visualize_tab.raster_figure

    def current_export_widget(self) -> QWidget:
        if self.currentIndex() == 1:
            return self.landscape_tab
        return self.stand_tab.current_export_widget()

    def current_figure(self) -> Figure:
        if self.currentIndex() == 1:
            return self.landscape_tab.figure
        return self.stand_tab.current_figure

    def current_unit_index(self) -> int:
        return self.stand_tab.current_unit_index()

    def current_unit(self) -> RxUnit:
        return self.stand_tab.current_unit()

    def current_raster_mode(self) -> int:
        return self.stand_tab.current_raster_mode()

    def showing_treatment_view(self) -> bool:
        return self.stand_tab.showing_treatment_view()

    def dbh_min(self) -> float:
        return self.stand_tab.dbh_min()

    def dbh_max(self) -> float:
        return self.stand_tab.dbh_max()

    def preview_dbh(self) -> float:
        return self.stand_tab.preview_dbh()

    def current_page(self) -> int:
        return self.stand_tab.current_page()

    def unit_system(self) -> UnitSystem:
        return self._state.unit_system

    def current_points(self) -> np.ndarray:
        return self.stand_tab.current_points()

    def current_raster_array(self) -> np.ndarray:
        return self.stand_tab.current_raster_array()

    def _handle_stand_state_change(self) -> None:
        self.landscape_tab.refresh()
        self._persistence.save_session(self._state)
