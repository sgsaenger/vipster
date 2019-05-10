cmake_minimum_required(VERSION 3.9)

project(Vipster VERSION 1.17)

option(TESTS "Build tests" OFF)
option(DESKTOP "Build QT-based desktop app" OFF)
option(BUILD_SHARED_LIBS "Build libvipster as shared library" ON)
include(CMakeDependentOption)
cmake_dependent_option(WEB "Build WASM-Module" OFF "NOT DESKTOP" OFF)
cmake_dependent_option(PYTHON "Python bindings" OFF "NOT WEB" OFF)
cmake_dependent_option(PYPI "To be set when built via setup.py" OFF "PYTHON" OFF)
cmake_dependent_option(EGG_INFO "Create infos for python package" ON "PYTHON" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(external)
include_directories(libvipster)
include_directories(python)
include_directories(vipster/qt)

# Library, needed for everything else
file(GLOB_RECURSE PY_SOURCES "libvipster/pyvipster.h" "libvipster/*.py.cpp")
file(GLOB_RECURSE LIB_SOURCES "libvipster/*.h" "libvipster/*.cpp" "external/tinyexpr.c")
list(REMOVE_ITEM LIB_SOURCES ${PY_SOURCES})
add_library(libvipster ${LIB_SOURCES})
configure_file("libvipster/version.h.in" "version.h")
set_target_properties(libvipster
    PROPERTIES
        FRAMEWORK TRUE
        OUTPUT_NAME "vipster"
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
if(UNIX AND NOT APPLE)
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

if(TESTS)
    # Tests, so far only for library
    enable_testing()
    file(GLOB TEST_SOURCES "tests/*.cpp")
    add_executable(test_lib ${TEST_SOURCES})
    target_link_libraries(test_lib libvipster)
    add_test(NAME test_lib COMMAND $<TARGET_FILE:test_lib>)
endif()

if(PYTHON)
    # Python-bindings
    find_package(pybind11 CONFIG)
    if(NOT pybind11_FOUND)
        add_subdirectory(external/pybind11)
    endif()
    pybind11_add_module(pyvipster ${PY_SOURCES})
    set_target_properties(pyvipster PROPERTIES OUTPUT_NAME vipster)
    target_link_libraries(pyvipster PRIVATE libvipster)
    if(PYPI OR (UNIX AND NOT APPLE))
        find_package(PythonInterp REQUIRED)
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
"from distutils import sysconfig as s;
print(s.get_python_lib(plat_specific=True,standard_lib=False))"
            OUTPUT_VARIABLE PYTHON_SITE
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        install(TARGETS pyvipster LIBRARY DESTINATION ${PYTHON_SITE}/vipster)
        install(FILES "util/package-init.py" DESTINATION ${PYTHON_SITE}/vipster RENAME "__init__.py")
        if(EGG_INFO AND NOT PYPI)
            install(CODE  "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" \"${PROJECT_SOURCE_DIR}/setup.py\" \"egg_info\")")
            install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vipster.egg-info" DESTINATION ${PYTHON_SITE})
        elseif(PYPI)
            configure_file("util/package-init.py" ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/__init__.py COPYONLY)
        endif()
    endif()
endif()

if(DESKTOP)
    # Qt-based desktop application

    # load and setup qt
    find_package(Qt5 CONFIG REQUIRED COMPONENTS Widgets Gui)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)

    # executable
    file(GLOB_RECURSE QT_SOURCES
         "vipster/qt/*.h"
         "vipster/qt/*.cpp"
         "vipster/qt/*.ui"
         "vipster/common/*.h"
         "vipster/common/*.cpp"
         "vipster/resources/vipster.qrc"
         )

    # icon-files:
    if(APPLE)
        set(MACOSX_BUNDLE_ICON_FILE vipster.icns)
        set(ICON_SOURCES "util/vipster.icns")
        set_source_files_properties(${ICON_SOURCES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    elseif(WIN32)
        set(ICON_SOURCES "util/win.rc")
    endif()
    add_executable(vipster WIN32 MACOSX_BUNDLE ${QT_SOURCES} ${ICON_SOURCES})
    target_link_libraries(vipster Qt5::Gui Qt5::Widgets libvipster)
    if(UNIX AND NOT APPLE)
        install(TARGETS vipster RUNTIME DESTINATION "bin")
        install(FILES "util/vipster.png" DESTINATION "share/icons/hicolor/128x128/apps")
        install(FILES "util/vipster.desktop" DESTINATION "share/applications")
        install(FILES "util/vipster.appdata.xml" DESTINATION "share/metainfo")
    endif()
endif()

if(WEB)
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
