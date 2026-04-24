# PyInstaller spec stub for RxGaming.
# This expects a Release build of rxgaming_core and a configured VCPKG_ROOT.

from __future__ import annotations

import os
from pathlib import Path

from PyInstaller.utils.hooks import Tree


ROOT = Path(SPECPATH).resolve().parent
PYTHON_DIR = ROOT / "python"
RESOURCES_DIR = ROOT / "resources"
ICON_PATH = ROOT / "icons" / "Icon.ico"
RELEASE_DIR = ROOT / "build" / "release" / "Release"
TRIPLET = os.environ.get("VCPKG_TARGET_TRIPLET", "x64-windows")
VCPKG_ROOT = Path(os.environ["VCPKG_ROOT"]) if "VCPKG_ROOT" in os.environ else None
VCPKG_BIN_DIR = VCPKG_ROOT / "installed" / TRIPLET / "bin" if VCPKG_ROOT else None


def existing_files(directory: Path, patterns: tuple[str, ...]) -> list[tuple[str, str]]:
    if not directory.exists():
        return []

    collected: list[tuple[str, str]] = []
    for pattern in patterns:
        for path in directory.glob(pattern):
            collected.append((str(path), "."))
    return collected


binaries = []
binaries += existing_files(PYTHON_DIR, ("rxgaming_core*.pyd",))
binaries += existing_files(RELEASE_DIR, ("*.dll",))

if VCPKG_BIN_DIR is not None:
    binaries += existing_files(VCPKG_BIN_DIR, ("*.dll",))

datas = []
if RESOURCES_DIR.exists():
    datas += Tree(str(RESOURCES_DIR), prefix="resources")


a = Analysis(
    ["python/__main__.py"],
    pathex=[str(ROOT), str(PYTHON_DIR)],
    binaries=binaries,
    datas=datas,
    hiddenimports=[
        "activity",
        "projectsettingsactivity",
        "gamingactivity",
        "gaming_tabs",
        "gaming_export",
        "widgets",
        "rxgaming_core",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="rxgaming",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    icon=str(ICON_PATH) if ICON_PATH.exists() else None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="rxgaming",
)
