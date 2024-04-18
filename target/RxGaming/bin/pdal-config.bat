@echo off

SET prefix=C:/OSGeo4W64
SET exec_prefix=C:/OSGeo4W64/bin
SET libdir=C:/OSGeo4W64/lib


IF "%1" == "--libs" echo -LC:/OSGeo4W64/lib -lpdalcpp & goto exit
IF "%1" == "--plugin-dir" echo C:/OSGeo4W64/bin & goto exit
IF "%1" == "--prefix" echo %prefix% & goto exit
IF "%1" == "--ldflags" echo -L%libdir% & goto exit
IF "%1" == "--defines" echo  & goto exit
IF "%1" == "--includes" echo -IC:/OSGeo4W64/include -IC:/OSGeo4W64/include -IC:/OSGeo4W64/include/libxml2 -IC:/OSGeo4W64/include -IC:/OSGeo4W64/include & goto exit
IF "%1" == "--cflags" echo /DWIN32 /D_WINDOWS /W3 & goto exit
IF "%1" == "--cxxflags" echo /DWIN32 /D_WINDOWS /W3 /GR /EHsc -std=c++11 & goto exit
IF "%1" == "--version" echo 1.8.0 & goto exit
IF "%1" == "--python-version" echo 3.7.0 & goto exit


echo Usage: pdal-config [OPTIONS]
echo Options:
echo    [--cflags]
echo    [--cxxflags]
echo    [--defines]
echo    [--includes]
echo    [--libs]
echo    [--plugin-dir]
echo    [--version]
echo    [--python-version]

:exit
