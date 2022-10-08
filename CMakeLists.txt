cmake_minimum_required(VERSION 3.24)

project(Vipster VERSION 1.20.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CMakeDependentOption)


#################
# CMake Options #
#################

# CMAKE_BUILD_TYPE - default to optimized Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel ..."
        FORCE)
endif()

# BUILD_TESTING - default to ON
option(BUILD_TESTING "Build tests" ON)

# BUILD_SHARED_LIBS - default to ON when not building for web
if(EMSCRIPTEN)
    # not supported with emscripten, hide warnings
    option(BUILD_SHARED_LIBS "Build shared library" OFF)
else()
    option(BUILD_SHARED_LIBS "Build shared library" ON)
endif()

# VIPSTER_PYLIB - add python bindings for library, default to OFF
find_package(Python3 3.7 QUIET COMPONENTS Interpreter Development)
cmake_dependent_option(VIPSTER_PYLIB "Standalone Python library" OFF ${Python3_FOUND} OFF)

# VIPSTER_DESKTOP - add Qt-GUI target, default to ON if Qt available
find_package(Qt6 CONFIG QUIET COMPONENTS Widgets Gui OpenGLWidgets)
if(NOT ${Qt6_FOUND})
    find_package(Qt5 5.10 CONFIG QUIET COMPONENTS Widgets Gui)
endif()
set(tmp ${VIPSTER_DESKTOP})
cmake_dependent_option(VIPSTER_DESKTOP "Build QT-based desktop app" ON "Qt6_FOUND OR Qt5_FOUND" OFF)
if(tmp AND NOT ${VIPSTER_DESKTOP})
    # ensure GUI is not silently disabled
    message(FATAL_ERROR "Requested to build the GUI, but no Qt installation was detected")
endif()

# VIPSTER_PYWIDGET - optional python widget for Qt GUI, default to ON if possible
cmake_dependent_option(VIPSTER_PYWIDGET "Python shell in desktop app" ${Python3_FOUND} "VIPSTER_DESKTOP" OFF)

# VIPSTER_LAMMPS - optional widget for interactive simulations, default to OFF
cmake_dependent_option(VIPSTER_LAMMPS "Interactive LAMMPS in desktop app" OFF "VIPSTER_DESKTOP" OFF)

# VIPSTER_WEB - minimal render widget for embedding in browser applications
cmake_dependent_option(VIPSTER_WEB "Build JS+WebGL GUI" ON "EMSCRIPTEN" OFF)

#################
# Build targets #
#################

# dependencies
include(util/dependencies.cmake)

# libvipster, core component that everything else depends on
add_subdirectory(vipster)

# Python bindings
if(VIPSTER_PYLIB OR VIPSTER_PYWIDGET)
    add_subdirectory(python)
endif()

# Full-featured Qt-based GUI
if(VIPSTER_DESKTOP)
    add_subdirectory(gui/qt)
endif()

# standalone browser widget
if(VIPSTER_WEB)
    add_subdirectory(gui/web)
endif()

# tests
if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
