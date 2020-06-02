cmake_minimum_required(VERSION 3.10)

project(Vipster VERSION 1.19)

# default to shared libraries
option(BUILD_SHARED_LIBS "Build shared library" ON)
# enable main GUI App if Qt is found
find_package(Qt5 5.10 CONFIG QUIET COMPONENTS Widgets Gui)
option(VIPSTER_DESKTOP "Build QT-based desktop app" ${Qt5_FOUND})
# automatic dependency solving
option(VIPSTER_DOWNLOAD_DEPENDENCIES "Download required dependencies if not present" ON)
mark_as_advanced(VIPSTER_DOWNLOAD_DEPENDENCIES)
include(CMakeDependentOption)
# Minimal standalone .wasm interface
cmake_dependent_option(VIPSTER_WEB "Build WASM-Module" OFF "NOT VIPSTER_DESKTOP" OFF)
# embedded LAMMPS widget
cmake_dependent_option(VIPSTER_LAMMPS "Interactive LAMMPS in desktop app" OFF "VIPSTER_DESKTOP" OFF)
# Python widget and standalone library
find_package(PythonInterp 3.6 QUIET)
cmake_dependent_option(VIPSTER_PYWIDGET "Python shell in desktop app" ${PythonInterp_FOUND} "VIPSTER_DESKTOP" OFF)
cmake_dependent_option(VIPSTER_PYLIB "Python library" ${PythonInterp_FOUND} "NOT VIPSTER_WEB" OFF)
cmake_dependent_option(VIPSTER_PYPI "To be set when built via setup.py" OFF "VIPSTER_PYLIB" OFF)
cmake_dependent_option(VIPSTER_EGG_INFO "Create infos for python package" ON "VIPSTER_PYLIB" OFF)
mark_as_advanced(VIPSTER_PYPI VIPSTER_EGG_INFO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Library, needed for everything else
file(GLOB_RECURSE PY_SOURCES "vipster/*.py.h" "vipster/*.py.cpp")
file(GLOB_RECURSE LIB_SOURCES "vipster/*.h" "vipster/*.cpp" "external/tinyexpr.c")
list(REMOVE_ITEM LIB_SOURCES ${PY_SOURCES})
list(REMOVE_ITEM PY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/vipster/main.py.cpp")
if(VIPSTER_WEB)
    # not needed here, so we can circumvent missing <filesystem>
    file(GLOB CONF_SOURCES "vipster/configfile.*")
    list(REMOVE_ITEM LIB_SOURCES ${CONF_SOURCES})
endif()
add_library(libvipster ${LIB_SOURCES})
configure_file("vipster/version.h.in" "version.h")
set_target_properties(libvipster
    PROPERTIES
        FRAMEWORK TRUE
        OUTPUT_NAME "vipster"
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
if(UNIX)
    # allow dynamic plugin loading
    target_link_libraries(libvipster PUBLIC dl)
    # on traditional unix-like platforms,
    # install files to appropriate locations
    if(BUILD_SHARED_LIBS AND NOT APPLE AND NOT VIPSTER_WEB)
        install(TARGETS libvipster
            LIBRARY
            DESTINATION "lib"
            FRAMEWORK
            DESTINATION "lib"
            )
        install(DIRECTORY "vipster/"
            DESTINATION "include/vipster"
            FILES_MATCHING PATTERN "*.h"
            )
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h DESTINATION "include/vipster")
    endif()
endif()

# function that intializes submodules if needed and requested
function(init_submodule NAME PATH)
    if(VIPSTER_DOWNLOAD_DEPENDENCIES)
        message("Could not find suitable ${NAME} installation, using included version at external/${PATH}")
        find_package(Git REQUIRED)
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --depth=1 external/${PATH}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    else()
        file(GLOB FILES ${CMAKE_CURRENT_SOURCE_DIR}/external/${PATH})
        list(LENGTH FILES NUM_FILES)
        if(NUM_FILES EQUAL 0)
            message(FATAL_ERROR
                "Could not find suitable ${NAME} installation and included version at external/${PATH} is not initialized.\n"
                "Please install ${NAME} or execute \"git submodule update --init ${PATH}\".")
        else()
            message("Could not find suitable ${NAME} installation, using included version at external/${PATH}")
        endif()
    endif()
endfunction()

# find suitable json-library
find_package(nlohmann_json 3.3.0 CONFIG QUIET)
if(NOT nlohmann_json_FOUND)
    init_submodule("nlohmann/json" "nlohmann_json")
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install OFF CACHE INTERNAL "")
    add_subdirectory(external/nlohmann_json EXCLUDE_FROM_ALL)
endif()

# find suitable fmt-library
find_package(fmt 6.1.2 CONFIG QUIET)
if(NOT fmt_FOUND)
    init_submodule("{fmt}" "fmt")
    add_subdirectory(external/fmt EXCLUDE_FROM_ALL)
endif()
target_link_libraries(libvipster PUBLIC nlohmann_json::nlohmann_json fmt::fmt-header-only)
target_include_directories(libvipster PUBLIC "external/tinyexpr")

if(BUILD_TESTING)
    # Tests, so far only for library
    file(GLOB TEST_SOURCES "tests/*.cpp")
    add_executable(test_lib ${TEST_SOURCES})

    # find suitable catch2-installation
    find_package(Catch2 2.4.1 CONFIG QUIET)
    if(NOT Catch2_FOUND)
        init_submodule("Catch2" "Catch2")
        set(CATCH_BUILD_TESTING OFF CACHE INTERNAL "")
        add_subdirectory(external/Catch2 EXCLUDE_FROM_ALL)
    endif()
    target_link_libraries(test_lib PRIVATE libvipster Catch2::Catch2)

    # enable testing via build system
    enable_testing()
    add_test(NAME test_lib COMMAND $<TARGET_FILE:test_lib>)
endif()

# create python binding objects
if(VIPSTER_PYLIB OR VIPSTER_PYWIDGET)
    find_package(PythonInterp 3.6 REQUIRED)
    find_package(PythonLibs 3.6 REQUIRED)
    find_package(pybind11 2.4.4 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        init_submodule("Pybind11" "pybind11")
        add_subdirectory(external/pybind11 EXCLUDE_FROM_ALL)
    endif()
    add_library(bindings OBJECT ${PY_SOURCES})
    target_link_libraries(bindings PUBLIC libvipster PRIVATE pybind11::module)
    if(WIN32)
        add_definitions(-DHAVE_SNPRINTF)
    endif()
endif()

if(VIPSTER_PYLIB)
    # standalone python bindings
    pybind11_add_module(pyvipster $<TARGET_OBJECTS:bindings> "vipster/main.py.cpp")
    set_target_properties(pyvipster PROPERTIES OUTPUT_NAME "vipster")
    target_link_libraries(pyvipster PRIVATE libvipster)
    if(VIPSTER_PYPI OR (UNIX AND NOT APPLE))
        execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
"from distutils import sysconfig as s;
print(s.get_python_lib(plat_specific=True,standard_lib=False))"
            OUTPUT_VARIABLE PYTHON_SITE
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        install(TARGETS pyvipster LIBRARY DESTINATION ${PYTHON_SITE})
        if(VIPSTER_EGG_INFO AND NOT VIPSTER_PYPI)
            install(CODE  "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" \"${PROJECT_SOURCE_DIR}/setup.py\" \"egg_info\")")
            install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vipster.egg-info" DESTINATION ${PYTHON_SITE})
        endif()
    endif()
endif()

if(VIPSTER_DESKTOP)
    # Qt-based desktop application

    # load and setup qt
    find_package(Qt5 5.10 CONFIG REQUIRED COMPONENTS Widgets Gui)
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)

    # find suitable cli11 installation
    find_package(CLI11 1.7.1 CONFIG QUIET)
    if(NOT CLI11_FOUND)
        init_submodule("CLI11" "CLI11")
        set(CLI11_TESTING OFF CACHE INTERNAL "")
        add_subdirectory(external/CLI11 EXCLUDE_FROM_ALL)
    endif()

    # get source files and split off optional sources
    file(GLOB_RECURSE QT_SOURCES
         "gui/qt/*.h"
         "gui/qt/*.cpp"
         "gui/qt/*.ui"
         "gui/common/*.h"
         "gui/common/*.cpp"
         "gui/resources/vipster.qrc"
         )
    file(GLOB_RECURSE QTPY_SOURCES "gui/*.py.*")
    list(REMOVE_ITEM QT_SOURCES ${QTPY_SOURCES})
    file(GLOB_RECURSE QTLMP_SOURCES "gui/*.lmp.*")
    list(REMOVE_ITEM QT_SOURCES ${QTLMP_SOURCES})

    # icon-files:
    if(APPLE)
        set(MACOSX_BUNDLE_ICON_FILE vipster.icns)
        set(ICON_SOURCES "util/vipster.icns")
        set_source_files_properties(${ICON_SOURCES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    elseif(WIN32)
        set(ICON_SOURCES "util/win.rc")
    endif()

    # executable
    add_executable(vipster WIN32 MACOSX_BUNDLE ${QT_SOURCES} ${ICON_SOURCES})

    # fix name-clash on import library when using mingw/gcc
    if(WIN32)
        set_target_properties(vipster PROPERTIES IMPORT_PREFIX "")
    endif()

    # link QT, CLI11 and vipster
    target_link_libraries(vipster PRIVATE Qt5::Gui Qt5::Widgets CLI11::CLI11 libvipster)
    target_include_directories(vipster PRIVATE "gui/common" "gui/qt")

    # python-widget:
    if(VIPSTER_PYWIDGET)
        add_definitions(-DUSE_PYTHON)
        target_sources(vipster PRIVATE ${QTPY_SOURCES} $<TARGET_OBJECTS:bindings>)
        target_link_libraries(vipster PRIVATE pybind11::embed)
    endif()

    # lammps-widget:
    if(VIPSTER_LAMMPS)
        set(USE_LAMMPS_BUILTIN ON)
        find_package(LAMMPS CONFIG QUIET)
        if(LAMMPS_FOUND)
            if(${LAMMPS_VERSION} GREATER_EQUAL 20200505)
                get_target_property(LMP_INC_DIR LAMMPS::lammps INTERFACE_INCLUDE_DIRECTORIES)
                if(EXISTS ${LMP_INC_DIR}/lammps/exceptions.h)
                    message("Found LAMMPS (suitable version \"${LAMMPS_VERSION}\" at ${LMP_INC_DIR})")
                    set(USE_LAMMPS_BUILTIN OFF)
                    # link against imported target
                    target_link_libraries(vipster PRIVATE LAMMPS::lammps)
                else()
                    message("Found incompatible LAMMPS (suitable version \"${LAMMPS_VERSION}\" at ${LMP_INC_DIR} has exceptions disabled)")
                endif()
            else()
                message("Found incompatible LAMMPS (version  \"${LAMMPS_VERSION}\" too low)")
            endif()
        endif()
        if(USE_LAMMPS_BUILTIN)
            init_submodule("LAMMPS" "lammps")
            # set LAMMPS build options as needed and hide them
            set(BUILD_TOOLS OFF CACHE INTERNAL "")
            set(BUILD_DOC OFF CACHE INTERNAL "")
            set(LAMMPS_EXCEPTIONS ON CACHE INTERNAL "")
            # enable simple and regularly used packages
            set(PKG_MOLECULE ON CACHE BOOL "")
            set(PKG_MANYBODY ON CACHE BOOL "")
            set(PKG_RIGID ON CACHE BOOL "")
            set(PKG_KSPACE ON CACHE BOOL "")
            set(PKG_MISC ON CACHE BOOL "")
            set(PKG_USER-MISC ON CACHE BOOL "")
            # include lammps as static library
            set(IS_SHARED BUILD_SHARED_LIBS)
            set(BUILD_SHARED_LIBS OFF)
            add_subdirectory(external/lammps/cmake EXCLUDE_FROM_ALL)
            set(BUILD_SHARED_LIBS IS_SHARED)
            # link against own built target (LAMMPS::lammps may be mis-configured by find_package)
            target_link_libraries(vipster PRIVATE lammps)
        endif()
        # enable widget in vipster
        add_definitions(-DUSE_LAMMPS)
        target_sources(vipster PRIVATE ${QTLMP_SOURCES})
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
    file(GLOB WEB_RESOURCES "gui/resources/*.frag" "gui/resources/*.vert")
    foreach(file ${WEB_RESOURCES})
        get_filename_component(filename ${file} NAME)
        set(EMCC_EMBED_FLAGS "${EMCC_EMBED_FLAGS}  --embed-file ${file}@${filename}")
    endforeach(file)
    set(CMAKE_EXE_LINKER_FLAGS "${EMCC_LINKER_FLAGS} ${EMCC_EMBED_FLAGS}")

    file(GLOB WEB_SOURCES
        "gh-pages/emscripten/index.html"
        "gh-pages/emscripten/vipster_setup.js"
        "gui/common/*.cpp"
        "gui/web/main.cpp")
    add_executable(webvipster ${WEB_SOURCES})
    set_target_properties(webvipster PROPERTIES OUTPUT_NAME vipster)
    target_link_libraries(webvipster PRIVATE libvipster)
    target_include_directories(webvipster PRIVATE "gui/common")
endif()
