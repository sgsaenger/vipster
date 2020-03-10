cmake_minimum_required(VERSION 3.9)

project(Vipster VERSION 1.18)

option(VIPSTER_TESTS "Build tests" OFF)
option(VIPSTER_DESKTOP "Build QT-based desktop app" OFF)
option(BUILD_SHARED_LIBS "Build libvipster as shared library" ON)
include(CMakeDependentOption)
cmake_dependent_option(VIPSTER_WEB "Build WASM-Module" OFF "NOT VIPSTER_DESKTOP" OFF)
cmake_dependent_option(VIPSTER_PYWIDGET "Python shell in desktop app" ON "VIPSTER_DESKTOP" OFF)
cmake_dependent_option(VIPSTER_PYLIB "Python bindings" OFF "NOT VIPSTER_WEB" OFF)
cmake_dependent_option(VIPSTER_PYPI "To be set when built via setup.py" OFF "VIPSTER_PYLIB" OFF)
cmake_dependent_option(VIPSTER_EGG_INFO "Create infos for python package" ON "VIPSTER_PYLIB" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(external/tinyexpr)
include_directories(libvipster)
include_directories(python)
include_directories(vipster/qt)

# Library, needed for everything else
file(GLOB_RECURSE PY_SOURCES "libvipster/pyvipster.h" "libvipster/*.py.h" "libvipster/*.py.cpp")
file(GLOB_RECURSE LIB_SOURCES "libvipster/*.h" "libvipster/*.cpp" "external/tinyexpr.c")
list(REMOVE_ITEM LIB_SOURCES ${PY_SOURCES})
list(REMOVE_ITEM PY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/libvipster/main.py.cpp")
if(VIPSTER_WEB)
    # not needed here, so we can circumvent missing <filesystem>
    file(GLOB CONF_SOURCES "libvipster/configfile.*")
    list(REMOVE_ITEM LIB_SOURCES ${CONF_SOURCES})
endif()
add_library(libvipster ${LIB_SOURCES})
configure_file("libvipster/version.h.in" "version.h")
set_target_properties(libvipster
    PROPERTIES
        FRAMEWORK TRUE
        OUTPUT_NAME "vipster"
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
if(UNIX)
    # allow dynamic plugin loading
    target_link_libraries(libvipster dl)
    # on traditional unix-like platforms,
    # install files to appropriate locations
    if(NOT APPLE AND NOT VIPSTER_WEB)
        install(TARGETS libvipster
            LIBRARY
            DESTINATION "lib"
            )
        install(DIRECTORY "libvipster/"
            DESTINATION "include/vipster"
            FILES_MATCHING PATTERN "*.h"
            )
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h DESTINATION "include/vipster")
    endif()
endif()
# find suitable json-library
find_package(nlohmann_json 3.3.0 CONFIG)
if(NOT nlohmann_json_FOUND)
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    add_subdirectory(external/nlohmann_json)
endif()
target_link_libraries(libvipster nlohmann_json::nlohmann_json)

if(VIPSTER_TESTS)
    # Tests, so far only for library
    file(GLOB TEST_SOURCES "tests/*.cpp")
    add_executable(test_lib ${TEST_SOURCES})

    # find suitable catch2-installation
    find_package(Catch2 2.4.1 CONFIG)
    if(NOT Catch2_FOUND)
        add_subdirectory(external/Catch2)
    endif()
    target_link_libraries(test_lib libvipster Catch2::Catch2)

    # enable testing via build system
    enable_testing()
    add_test(NAME test_lib COMMAND $<TARGET_FILE:test_lib>)
endif()

# create python binding objects
if(VIPSTER_PYLIB OR VIPSTER_PYWIDGET)
    find_package(PythonInterp 3.6 REQUIRED)
    find_package(PythonLibs 3.6 REQUIRED)
    find_package(pybind11 2.4.4 CONFIG)
    if(NOT pybind11_FOUND)
        add_subdirectory(external/pybind11)
    endif()
    add_library(pyvipster OBJECT ${PY_SOURCES})
    target_link_libraries(pyvipster PRIVATE pybind11::module)
endif()

if(VIPSTER_PYLIB)
    # standalone python bindings
    pybind11_add_module(bindings $<TARGET_OBJECTS:pyvipster> "libvipster/main.py.cpp")
    set_target_properties(bindings PROPERTIES OUTPUT_NAME "_vipster")
    target_link_libraries(bindings PRIVATE libvipster)
    if(VIPSTER_PYPI OR (UNIX AND NOT APPLE))
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
"from distutils import sysconfig as s;
print(s.get_python_lib(plat_specific=True,standard_lib=False))"
            OUTPUT_VARIABLE PYTHON_SITE
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        install(TARGETS pyvipster LIBRARY DESTINATION ${PYTHON_SITE}/vipster)
        install(FILES "util/package-init.py" DESTINATION ${PYTHON_SITE}/vipster RENAME "__init__.py")
        if(VIPSTER_EGG_INFO AND NOT VIPSTER_PYPI)
            install(CODE  "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" \"${PROJECT_SOURCE_DIR}/setup.py\" \"egg_info\")")
            install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vipster.egg-info" DESTINATION ${PYTHON_SITE})
        elseif(VIPSTER_PYPI)
            configure_file("util/package-init.py" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/__init__.py COPYONLY)
        endif()
    endif()
endif()

if(VIPSTER_DESKTOP)
    # Qt-based desktop application

    # load and setup qt
    find_package(Qt5 5.10 CONFIG REQUIRED COMPONENTS Widgets Gui)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)

    # find suitable cli11 installation
    find_package(CLI11 1.7.1 CONFIG)
    if(NOT CLI11_FOUND)
        set(CLI11_TESTING OFF CACHE INTERNAL "")
        add_subdirectory(external/CLI11)
    endif()

    # executable
    file(GLOB_RECURSE QT_SOURCES
         "vipster/qt/*.h"
         "vipster/qt/*.cpp"
         "vipster/qt/*.ui"
         "vipster/common/*.h"
         "vipster/common/*.cpp"
         "vipster/resources/vipster.qrc"
         )

    # python-widget:
    if(NOT VIPSTER_PYWIDGET)
        file(GLOB_RECURSE QTPY_SOURCES "vipster/*.py.*")
        list(REMOVE_ITEM QT_SOURCES ${QTPY_SOURCES})
    else()
        add_definitions(-DUSE_PYTHON)
        list(APPEND QT_SOURCES $<TARGET_OBJECTS:pyvipster>)
    endif()

    # icon-files:
    if(APPLE)
        set(MACOSX_BUNDLE_ICON_FILE vipster.icns)
        set(ICON_SOURCES "util/vipster.icns")
        set_source_files_properties(${ICON_SOURCES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    elseif(WIN32)
        set(ICON_SOURCES "util/win.rc")
    endif()
    add_executable(vipster WIN32 MACOSX_BUNDLE ${QT_SOURCES} ${ICON_SOURCES})
    # fix name-clash on import lib, assuming gcc-conventions, needs changes for clang
    if(WIN32)
        set_target_properties(vipster PROPERTIES IMPORT_PREFIX "")
    endif()
    # set link targets
    target_link_libraries(vipster Qt5::Gui Qt5::Widgets CLI11::CLI11 libvipster)
    if(VIPSTER_PYWIDGET)
        target_link_libraries(vipster pybind11::embed)
    endif()
    # on traditional unix-like platforms,
    # install files to appropriate locations
    if(UNIX AND NOT APPLE)
        install(TARGETS vipster RUNTIME DESTINATION "bin")
        install(FILES "util/vipster.png" DESTINATION "share/icons/hicolor/128x128/apps")
        install(FILES "util/vipster.desktop" DESTINATION "share/applications")
        install(FILES "util/vipster.appdata.xml" DESTINATION "share/metainfo")
    endif()
endif()

if(VIPSTER_WEB)
    # Emscripten-based web application
    set(EMCC_LINKER_FLAGS "--bind -s USE_WEBGL2=1 -s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0")
    file(GLOB WEB_RESOURCES "vipster/resources/*.frag" "vipster/resources/*.vert")
    foreach(file ${WEB_RESOURCES})
        get_filename_component(filename ${file} NAME)
        set(EMCC_EMBED_FLAGS "${EMCC_EMBED_FLAGS}  --embed-file ${file}@${filename}")
    endforeach(file)
    set(CMAKE_EXE_LINKER_FLAGS "${EMCC_LINKER_FLAGS} ${EMCC_EMBED_FLAGS}")

    file(GLOB WEB_SOURCES
        "gh-pages/emscripten/index.html"
        "gh-pages/emscripten/vipster_setup.js"
        "vipster/common/*.cpp"
        "vipster/web/main.cpp")
    add_executable(webvipster ${WEB_SOURCES})
    set_target_properties(webvipster PROPERTIES OUTPUT_NAME vipster)
    target_link_libraries(webvipster libvipster)
endif()
