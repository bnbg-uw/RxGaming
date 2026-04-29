from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.figure import Figure
import numpy as np
from PySide6.QtWidgets import QFileDialog, QMessageBox, QWidget

from gaming_ui import GamingTabs
from gaming_ui.units import UnitSystem, array_to_display, display_name_for


_BASIN_PALETTE_SIZE = 64
_CLUMP_COLORS = ("white", "#7bc043", "#fdf498", "#f37736", "#ee4035")
_CLUMP_BOUNDARIES = (-0.5, 0.5, 1.5, 4.5, 9.5, 99)
_CLUMP_TICKS = [0, 1, 3, 7, 55]
_CLUMP_TICK_LABELS = ["0", "1", "2-4", "4-9", "10+"]


@dataclass(frozen=True)
class _GeoRasterOption:
    label: str
    output_stem: str
    render_mode: str
    data_method: str
    hillshade_method: str | None = None


def export_current_view(tabs: GamingTabs, parent: QWidget) -> None:
    file_path = QFileDialog.getSaveFileName(parent, "Export current view", "", "PNG files (*.png)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() != ".png":
        output_path = output_path.with_suffix(".png")

    tabs.current_export_widget().grab().save(str(output_path), "PNG")


def export_georeferenced_raster(tabs: GamingTabs, parent: QWidget) -> None:
    try:
        unit = tabs.current_unit()
    except (AttributeError, IndexError):
        _show_warning(parent, "No unit is available to export.")
        return

    option = _current_georeferenced_raster_option(tabs)
    default_name = f"{_safe_output_name(getattr(unit, 'name', 'unit'))}_{option.output_stem}.tif"
    file_path = QFileDialog.getSaveFileName(parent, "Export georeferenced raster image", default_name, "TIFF files (*.tif *.tiff)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() not in {".tif", ".tiff"}:
        output_path = output_path.with_suffix(".tif")

    export_method = getattr(unit, "export_rendered_geotiff", None)
    if export_method is None:
        _show_warning(parent, f"Native export support for {option.label} is not available.")
        return

    render = _render_georeferenced_raster(unit, option, tabs.unit_system())
    try:
        export_method(
            str(output_path),
            render.image,
            render.map_left_px,
            render.map_top_px,
            render.map_width_px,
            render.map_height_px,
        )
    except Exception as exc:
        _show_warning(parent, f"Could not export georeferenced raster image.\n\n{exc}")


def export_chm_raster(tabs: GamingTabs, parent: QWidget) -> None:
    _export_native_raster(
        tabs,
        parent,
        title="Export CHM",
        current_label="Current CHM",
        treated_label="Treated CHM",
        output_stem="current_chm",
        treated_output_stem="treated_chm",
        method_name="write_chm_raster",
    )


def export_basins_raster(tabs: GamingTabs, parent: QWidget) -> None:
    _export_native_raster(
        tabs,
        parent,
        title="Export Basins",
        current_label="Current Basin Map",
        treated_label="Treated Basin Map",
        output_stem="current_basin",
        treated_output_stem="treated_basin",
        method_name="write_basin_raster",
    )


def export_clumpmap_raster(tabs: GamingTabs, parent: QWidget) -> None:
    _export_native_raster(
        tabs,
        parent,
        title="Export Clumpmap",
        current_label="Current Clump Map",
        treated_label="Treated Clump Map",
        output_stem="current_clump",
        treated_output_stem="treated_clump",
        method_name="write_clumpmap_raster",
    )


def export_features(tabs: GamingTabs, parent: QWidget) -> None:
    try:
        unit = tabs.current_unit()
    except (AttributeError, IndexError):
        _show_warning(parent, "No unit is available to export.")
        return

    default_name = f"{_safe_output_name(getattr(unit, 'name', 'unit'))}_points.shp"
    file_path = QFileDialog.getSaveFileName(parent, "Export point data", default_name, "ESRI Shapefile (*.shp)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() != ".shp":
        output_path = output_path.with_suffix(".shp")

    export_method = getattr(unit, "write_tao_shapefile", None)
    if export_method is None:
        _show_warning(parent, "Native export support for point data is not available.")
        return

    try:
        export_method(str(output_path), tabs.showing_treatment_view())
    except Exception as exc:
        _show_warning(parent, f"Could not export point data.\n\n{exc}")


def export_treelist(tabs: GamingTabs, parent: QWidget) -> None:
    _export_points_csv(tabs, parent, "Export treelist", output_stem="treelist")


@dataclass(frozen=True)
class _RenderedGeoRaster:
    image: np.ndarray
    map_left_px: int
    map_top_px: int
    map_width_px: int
    map_height_px: int


def _current_georeferenced_raster_option(tabs: GamingTabs) -> _GeoRasterOption:
    raster_mode = tabs.current_raster_mode()
    show_treatment = tabs.showing_treatment_view()
    if raster_mode == 0:
        return _GeoRasterOption(
            "Treated CHM" if show_treatment else "Current CHM",
            "treated_chm_map" if show_treatment else "current_chm",
            "chm",
            "get_treat_chm" if show_treatment else "get_chm",
            "get_treat_hillshade" if show_treatment else "get_hillshade",
        )
    if raster_mode == 1:
        return _GeoRasterOption(
            "Treated Basin Map" if show_treatment else "Current Basin Map",
            "treated_basin_map" if show_treatment else "current_basin",
            "basin",
            "get_treat_basin" if show_treatment else "get_basin",
            "get_treat_hillshade" if show_treatment else "get_hillshade",
        )
    if raster_mode == 2:
        return _GeoRasterOption(
            "Treated Clump Map" if show_treatment else "Current Clump Map",
            "treated_clump_map" if show_treatment else "current_clump",
            "clump",
            "get_treat_clump_map" if show_treatment else "get_clump_map",
            "get_treat_hillshade" if show_treatment else "get_hillshade",
        )
    return _GeoRasterOption(
        "Treated Hillshade" if show_treatment else "Current Hillshade",
        "treated_hillshade_map" if show_treatment else "current_hillshade",
        "hillshade",
        "get_treat_hillshade" if show_treatment else "get_hillshade",
    )


def _render_georeferenced_raster(unit: object, option: _GeoRasterOption, unit_system: UnitSystem) -> _RenderedGeoRaster:
    data = np.asarray(getattr(unit, option.data_method)())
    hillshade = np.asarray(getattr(unit, option.hillshade_method)()) if option.hillshade_method is not None else np.zeros((0, 0))

    figure = Figure(figsize=(10.0, 7.5), dpi=300)
    canvas = FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    figure.subplots_adjust(left=0.03, right=0.84, top=0.94, bottom=0.04)
    colorbar_axes = figure.add_axes([0.86, 0.04, 0.03, 0.90])
    colorbar_axes.set_visible(False)

    image_artist = None
    colorbar_label = None
    colorbar_ticks = None
    colorbar_ticklabels = None

    if option.render_mode == "chm":
        display_data = array_to_display("canopy_height", data, unit_system)
        image_artist = axes.imshow(display_data, cmap="coolwarm", vmin=0)
        if hillshade.size:
            axes.imshow(hillshade, cmap="Greys", alpha=0.5)
        _set_dynamic_limits(image_artist, display_data, lower_bound=0.0)
        colorbar_label = display_name_for("canopy_height", unit_system)
    elif option.render_mode == "basin":
        basin_display = _categorical_basin_display(data, option.data_method.startswith("get_treat_"))
        image_artist = axes.imshow(
            basin_display,
            cmap=_build_basin_colormap(),
            interpolation="nearest",
            vmin=0,
            vmax=_BASIN_PALETTE_SIZE - 1,
        )
        if hillshade.size:
            axes.imshow(hillshade, cmap="Greys", alpha=0.5)
    elif option.render_mode == "clump":
        cmap = LinearSegmentedColormap.from_list("Clump Colors", _CLUMP_COLORS, 5)
        norm = BoundaryNorm(_CLUMP_BOUNDARIES, len(_CLUMP_BOUNDARIES))
        image_artist = axes.imshow(data, cmap=cmap, norm=norm)
        if hillshade.size:
            axes.imshow(hillshade, cmap="Greys", alpha=0.5)
        colorbar_label = "Clump Map (Clump bins)"
        colorbar_ticks = _CLUMP_TICKS
        colorbar_ticklabels = _CLUMP_TICK_LABELS
    else:
        image_artist = axes.imshow(data, cmap="Greys")
        _set_dynamic_limits(image_artist, data)
        colorbar_label = "Hillshade"

    axes.set_title(f"{getattr(unit, 'name', 'Unit')} {option.label}")
    axes.set_xticks([])
    axes.set_yticks([])

    if image_artist is not None and colorbar_label is not None:
        colorbar_axes.set_visible(True)
        colorbar = figure.colorbar(image_artist, cax=colorbar_axes)
        colorbar.set_label(colorbar_label)
        if colorbar_ticks is not None:
            colorbar.set_ticks(colorbar_ticks)
        if colorbar_ticklabels is not None:
            colorbar.set_ticklabels(colorbar_ticklabels)

    canvas.draw()
    width_px, height_px = canvas.get_width_height()
    image = np.asarray(canvas.buffer_rgba()).copy()
    bbox = axes.get_window_extent(renderer=canvas.get_renderer())
    map_left_px = max(0, int(round(bbox.x0)))
    map_top_px = max(0, int(round(height_px - bbox.y1)))
    map_width_px = max(1, int(round(bbox.width)))
    map_height_px = max(1, int(round(bbox.height)))

    return _RenderedGeoRaster(
        image=image,
        map_left_px=map_left_px,
        map_top_px=map_top_px,
        map_width_px=map_width_px,
        map_height_px=map_height_px,
    )


def _export_points_csv(tabs: GamingTabs, parent: QWidget, title: str, output_stem: str = "points") -> None:
    default_name = ""
    try:
        unit = tabs.current_unit()
    except (AttributeError, IndexError):
        unit = None
    if unit is not None:
        default_name = f"{_safe_output_name(getattr(unit, 'name', 'unit'))}_{output_stem}.csv"

    file_path = QFileDialog.getSaveFileName(parent, title, default_name, "CSV files (*.csv)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() != ".csv":
        output_path = output_path.with_suffix(".csv")

    points = np.asarray(tabs.current_points())
    if points.ndim != 2 or points.shape[1] == 0:
        _show_warning(parent, "No point data is available to export.")
        return

    headers = _point_csv_headers(points.shape[1])
    with output_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(headers)
        writer.writerows(points.tolist())


def _point_csv_headers(column_count: int) -> list[str]:
    if column_count == 2:
        return ["x", "y"]
    if column_count == 4:
        return ["x", "y", "height", "radius"]
    if column_count == 5:
        return ["x", "y", "height", "radius", "dbh"]
    if column_count == 6:
        return ["x", "y", "area", "height", "radius", "dbh"]
    return [f"column_{index + 1}" for index in range(column_count)]


def _export_native_raster(
    tabs: GamingTabs,
    parent: QWidget,
    *,
    title: str,
    current_label: str,
    treated_label: str,
    output_stem: str,
    treated_output_stem: str,
    method_name: str,
) -> None:
    try:
        unit = tabs.current_unit()
    except (AttributeError, IndexError):
        _show_warning(parent, "No unit is available to export.")
        return

    show_treatment = tabs.showing_treatment_view()
    label = treated_label if show_treatment else current_label
    selected_output_stem = treated_output_stem if show_treatment else output_stem
    default_name = f"{_safe_output_name(getattr(unit, 'name', 'unit'))}_{selected_output_stem}.tif"
    file_path = QFileDialog.getSaveFileName(parent, title, default_name, "TIFF files (*.tif *.tiff)")[0]
    if not file_path:
        return

    output_path = Path(file_path)
    if output_path.suffix.lower() not in {".tif", ".tiff"}:
        output_path = output_path.with_suffix(".tif")

    export_method = getattr(unit, method_name, None)
    if export_method is None:
        _show_warning(parent, f"Native export support for {label} is not available.")
        return

    try:
        export_method(str(output_path), show_treatment)
    except Exception as exc:
        _show_warning(parent, f"Could not export {label.lower()}.\n\n{exc}")


def _safe_output_name(name: str) -> str:
    sanitized = "".join(character if character.isalnum() or character in {"-", "_"} else "_" for character in name.strip())
    return sanitized or "unit"


def _build_basin_colormap() -> object:
    from matplotlib.colors import ListedColormap, hsv_to_rgb

    hues = (np.arange(_BASIN_PALETTE_SIZE, dtype=float) * 0.6180339887498949) % 1.0
    saturations = np.where(np.arange(_BASIN_PALETTE_SIZE) % 2 == 0, 0.70, 0.55)
    values = np.where(np.arange(_BASIN_PALETTE_SIZE) % 4 < 2, 0.92, 0.78)
    hsv = np.column_stack((hues, saturations, values))
    colors = hsv_to_rgb(hsv)
    colors[0] = (0.65, 0.65, 0.65)
    cmap = ListedColormap(colors, name="Basin Colors")
    cmap.set_bad((1.0, 1.0, 1.0, 1.0))
    return cmap


def _categorical_basin_display(data: np.ndarray, show_treatment: bool) -> np.ma.MaskedArray:
    array = np.asarray(data)
    if array.size == 0:
        return np.ma.masked_array(array, mask=np.ones_like(array, dtype=bool))

    if not np.issubdtype(array.dtype, np.integer):
        integer_array = array.astype(np.int64, copy=False)
    else:
        integer_array = array

    dtype_sentinel = np.iinfo(integer_array.dtype).min
    int32_sentinel = np.int64(np.iinfo(np.int32).min)
    mask = (integer_array == dtype_sentinel) | (integer_array == int32_sentinel)
    hashed = np.zeros(integer_array.shape, dtype=np.uint8)
    removed_mask = (~mask) & (integer_array == 1) if show_treatment else np.zeros_like(mask, dtype=bool)
    valid_tree_mask = (~mask) & (~removed_mask)
    if np.any(valid_tree_mask):
        ids = integer_array[valid_tree_mask].astype(np.uint64, copy=False)
        ids ^= ids >> np.uint64(33)
        ids *= np.uint64(0xff51afd7ed558ccd)
        ids ^= ids >> np.uint64(33)
        ids *= np.uint64(0xc4ceb9fe1a85ec53)
        ids ^= ids >> np.uint64(33)
        hashed[valid_tree_mask] = (np.mod(ids, _BASIN_PALETTE_SIZE - 1) + 1).astype(np.uint8, copy=False)
    hashed[removed_mask] = 0
    return np.ma.masked_array(hashed, mask=mask)


def _set_dynamic_limits(image_artist: object, data: np.ndarray, lower_bound: float | None = None) -> None:
    min_value, max_value = _array_bounds(data)
    if lower_bound is not None:
        min_value = lower_bound
    if min_value is None or max_value is None:
        min_value, max_value = (0.0, 1.0)
    if max_value <= min_value:
        max_value = min_value + 1.0
    image_artist.set_clim(float(min_value), float(max_value))


def _array_bounds(data: np.ndarray) -> tuple[float | None, float | None]:
    array = np.asarray(data)
    if array.size == 0:
        return (None, None)
    if np.issubdtype(array.dtype, np.integer):
        sentinel = np.iinfo(array.dtype).min
        array = array[array != sentinel]
    else:
        array = array[np.isfinite(array)]
    if array.size == 0:
        return (None, None)
    return (float(np.min(array)), float(np.max(array)))


def _show_warning(parent: QWidget, text: str) -> None:
    msg = QMessageBox(parent)
    msg.setIcon(QMessageBox.Icon.Warning)
    msg.setText(text)
    msg.setWindowTitle("Export unavailable")
    msg.setStandardButtons(QMessageBox.StandardButton.Ok)
    msg.exec()
