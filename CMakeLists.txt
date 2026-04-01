cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0148 NEW)
cmake_policy(SET CMP0167 NEW)

include_guard(GLOBAL)
project(rxgaming)

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

# Submodule dependencies
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

target_include_directories(rxgaming_core PRIVATE src/)

target_link_libraries(rxgaming_core PRIVATE
    HDF5::HDF5
    nlohmann_json::nlohmann_json
    OpenMP::OpenMP_CXX
    RxTools
)

# Copy built module to python/ folder
add_custom_command(TARGET rxgaming_core POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:rxgaming_core>
        ${CMAKE_SOURCE_DIR}/python/
)

target_compile_features(rxgaming_core PRIVATE cxx_std_20)
target_compile_options(rxgaming_core PRIVATE
    $<$<CXX_COMPILER_ID:MSVC>:/W3 /WX>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Werror>
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)