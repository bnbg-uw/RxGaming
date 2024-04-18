# -*- mode: python ; coding: utf-8 -*-


a = Analysis(
    ['F:\\RxGaming\\src\\main\\python\\rxgaming\\__main__.py'],
    pathex=[],
    binaries=[],
    datas=[],
    hiddenimports=['rasterio._base', 'rasterio._crs', 'rasterio._env', 'rasterio._err', 'rasterio._example', 'rasterio._features', 'rasterio._fill', 'rasterio._io', 'rasterio._shim', 'rasterio._warp', 'rasterio.compat', 'rasterio.control', 'rasterio.coords', 'rasterio.crs', 'rasterio.drivers', 'rasterio.dtypes', 'rasterio.enums', 'rasterio.env', 'rasterio.errors', 'rasterio.features', 'rasterio.fill', 'rasterio.io', 'rasterio.mask', 'rasterio.merge', 'rasterio.path', 'rasterio.plot', 'rasterio.profiles', 'rasterio.sample', 'rasterio.session', 'rasterio.shutil', 'rasterio.tools', 'rasterio.transform', 'rasterio.vrt', 'rasterio.warp', 'rasterio.windows', 'osgeo._gdal', 'osgeo.gdal', 'richdem.__init__', 'richdem.cli'],
    hookspath=['F:\\RxGaming\\.venv\\lib\\site-packages\\fbs\\freeze\\hooks'],
    hooksconfig={},
    runtime_hooks=['F:\\RxGaming\\target\\PyInstaller\\fbs_pyinstaller_hook.py'],
    excludes=[],
    noarchive=True,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [('v', None, 'OPTION')],
    exclude_binaries=True,
    name='RxGaming',
    debug=True,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    version='F:\\RxGaming\\target\\PyInstaller\\version_info.py',
    icon=['F:\\RxGaming\\src\\main\\icons\\Icon.ico'],
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='RxGaming',
)
