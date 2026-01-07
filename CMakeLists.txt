cmake_minimum_required(VERSION 3.12)

project(MicroUGens
    VERSION 1.0.0
    DESCRIPTION "Microsound UGens for SuperCollider"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Options
option(SUPERNOVA "Build plugins for supernova" OFF)
option(SC_PATH "Path to SuperCollider source" "")

# Find SuperCollider source
if(SC_PATH)
    set(SC_SOURCE_DIR ${SC_PATH})
elseif(DEFINED ENV{SC_PATH})
    set(SC_SOURCE_DIR $ENV{SC_PATH})
else()
    if(EXISTS "${CMAKE_SOURCE_DIR}/../supercollider")
        set(SC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/../supercollider")
    elseif(EXISTS "${CMAKE_SOURCE_DIR}/supercollider")
        set(SC_SOURCE_DIR "${CMAKE_SOURCE_DIR}/supercollider")
    elseif(EXISTS "$ENV{HOME}/supercollider")
        set(SC_SOURCE_DIR "$ENV{HOME}/supercollider")
    elseif(EXISTS "$ENV{HOME}/src/supercollider")
        set(SC_SOURCE_DIR "$ENV{HOME}/src/supercollider")
    else()
        message(FATAL_ERROR
            "SuperCollider source not found. Set SC_PATH:\n"
            "  cmake -DSC_PATH=/path/to/supercollider ..")
    endif()
endif()

message(STATUS "SuperCollider source: ${SC_SOURCE_DIR}")

if(NOT EXISTS "${SC_SOURCE_DIR}/include/plugin_interface/SC_PlugIn.h")
    message(FATAL_ERROR "SC_PlugIn.h not found in ${SC_SOURCE_DIR}/include/plugin_interface")
endif()

# Include directories
include_directories(
    ${SC_SOURCE_DIR}/include/plugin_interface
    ${SC_SOURCE_DIR}/include/common
    ${SC_SOURCE_DIR}/common
)

# Platform settings
if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS version")
    set(CMAKE_OSX_ARCHITECTURES "x86_64;arm64" CACHE STRING "Build architectures")
    set(PLUGIN_EXTENSION ".scx")
    add_compile_definitions(SC_DARWIN)
    set(INSTALL_PATH "$ENV{HOME}/Library/Application Support/SuperCollider/Extensions/MicroUGens")
elseif(WIN32)
    set(PLUGIN_EXTENSION ".scx")
    add_compile_definitions(SC_WIN32)
    set(INSTALL_PATH "$ENV{LOCALAPPDATA}/SuperCollider/Extensions/MicroUGens")
else()
    set(PLUGIN_EXTENSION ".so")
    add_compile_definitions(SC_LINUX)
    add_compile_options(-fPIC)
    set(INSTALL_PATH "$ENV{HOME}/.local/share/SuperCollider/Extensions/MicroUGens")
endif()

# Compiler flags
add_compile_options(-O3 -ffast-math)
if(NOT APPLE)
    add_compile_options(-march=native)
endif()

# --- SubGrain ---
add_library(SubGrain MODULE plugins/SubGrain/SubGrain.cpp)
set_target_properties(SubGrain PROPERTIES PREFIX "" SUFFIX ${PLUGIN_EXTENSION})

# --- Pulsar ---
add_library(Pulsar MODULE plugins/Pulsar/Pulsar.cpp)
set_target_properties(Pulsar PROPERTIES PREFIX "" SUFFIX ${PLUGIN_EXTENSION})

# Supernova builds (optional)
if(SUPERNOVA)
    add_library(SubGrain_supernova MODULE plugins/SubGrain/SubGrain.cpp)
    target_compile_definitions(SubGrain_supernova PRIVATE SUPERNOVA)
    set_target_properties(SubGrain_supernova PROPERTIES PREFIX "" SUFFIX ${PLUGIN_EXTENSION})

    add_library(Pulsar_supernova MODULE plugins/Pulsar/Pulsar.cpp)
    target_compile_definitions(Pulsar_supernova PRIVATE SUPERNOVA)
    set_target_properties(Pulsar_supernova PROPERTIES PREFIX "" SUFFIX ${PLUGIN_EXTENSION})
endif()

# Install
install(TARGETS SubGrain Pulsar LIBRARY DESTINATION ${INSTALL_PATH})
install(FILES plugins/SubGrain/SubGrain.sc plugins/Pulsar/Pulsar.sc
        DESTINATION ${INSTALL_PATH}/Classes)
install(FILES plugins/SubGrain/SubGrain.schelp plugins/Pulsar/Pulsar.schelp
        DESTINATION ${INSTALL_PATH}/HelpSource/Classes)

if(SUPERNOVA)
    install(TARGETS SubGrain_supernova Pulsar_supernova LIBRARY DESTINATION ${INSTALL_PATH})
endif()

# Summary
message(STATUS "")
message(STATUS "MicroUGens configuration:")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME}")
message(STATUS "  Supernova: ${SUPERNOVA}")
message(STATUS "  Install path: ${INSTALL_PATH}")
message(STATUS "")
message(STATUS "Build:")
message(STATUS "  mkdir build && cd build")
message(STATUS "  cmake -DSC_PATH=/path/to/supercollider ..")
message(STATUS "  make && make install")
message(STATUS "")
