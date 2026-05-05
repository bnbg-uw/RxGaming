cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0148 NEW)
cmake_policy(SET CMP0167 NEW)

include_guard(GLOBAL)
project(rxgaming)
enable_testing()

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

# Submodule dependencies
set(LAPISGIS_AUTO_INIT OFF CACHE BOOL "Call lapisGisInit before main" FORCE)
add_subdirectory(src/rxtools)

# vcpkg dependencies
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(nlohmann_json REQUIRED)

set(PYBIND11_FINDPYTHON ON)
find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

find_package(OpenMP REQUIRED)

# Your source files
file(GLOB_RECURSE SOURCES "src/*.cpp")
list(FILTER SOURCES EXCLUDE REGEX "src/(lapisgis|lico|processedfolder|rxtools)/.*")

pybind11_add_module(rxgaming_core ${SOURCES})

set(RXGAMING_CORE_CPP_SOURCES ${SOURCES})
list(FILTER RXGAMING_CORE_CPP_SOURCES EXCLUDE REGEX "src/main_lib\\.cpp$")

target_include_directories(rxgaming_core PRIVATE src/)

target_link_libraries(rxgaming_core PRIVATE
    HDF5::HDF5
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
    RxTools
)

add_executable(rxgaming_projectarea_smoke
    tools/projectarea_smoke.cpp
    ${RXGAMING_CORE_CPP_SOURCES}
)

target_include_directories(rxgaming_projectarea_smoke PRIVATE src/)

target_link_libraries(rxgaming_projectarea_smoke PRIVATE
    HDF5::HDF5
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
    Python::Python
    pybind11::headers
    RxTools
)

add_test(
    NAME rxgaming_projectarea_smoke_help
    COMMAND rxgaming_projectarea_smoke --help
)

# Copy built module to python/ folder
add_custom_command(TARGET rxgaming_core POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:rxgaming_core>
        ${CMAKE_SOURCE_DIR}/python/
)

set(PROJ_DB_PATH "" CACHE FILEPATH "Path to proj.db")
if(NOT PROJ_DB_PATH)
    set(_proj_db_candidates "")

    if(DEFINED ENV{VCPKG_ROOT} AND DEFINED VCPKG_TARGET_TRIPLET)
        list(APPEND _proj_db_candidates
            "$ENV{VCPKG_ROOT}/installed/${VCPKG_TARGET_TRIPLET}/share/proj/proj.db"
        )
    endif()

    foreach(_proj_db_candidate IN LISTS _proj_db_candidates)
        if(EXISTS "${_proj_db_candidate}")
            set(PROJ_DB_PATH "${_proj_db_candidate}" CACHE FILEPATH "Path to proj.db" FORCE)
            message(STATUS "Auto-detected proj.db at: ${PROJ_DB_PATH}")
            break()
        endif()
    endforeach()
endif()

if(EXISTS "${PROJ_DB_PATH}")
    add_custom_command(TARGET rxgaming_core POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_SOURCE_DIR}/resources
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${PROJ_DB_PATH}"
            "${CMAKE_SOURCE_DIR}/resources/proj.db"
        COMMENT "Copying proj.db to resources/"
    )
else()
    message(WARNING "proj.db not found at: ${PROJ_DB_PATH}\nSet -DPROJ_DB_PATH to copy proj.db into resources/.")
endif()

target_compile_features(rxgaming_core PRIVATE cxx_std_20)
target_compile_options(rxgaming_core PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Werror>
)

target_compile_features(rxgaming_projectarea_smoke PRIVATE cxx_std_20)
target_compile_options(rxgaming_projectarea_smoke PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Werror>
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
