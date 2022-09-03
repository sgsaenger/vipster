cmake_minimum_required(VERSION 3.14)

project(Vipster VERSION 1.19.1)

include(CMakeDependentOption)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel ..."
        FORCE)
endif()

# automatic dependency solving
option(VIPSTER_DOWNLOAD_DEPENDENCIES "Download required dependencies if not present" ON)

# detect compilation with emscripten
if(EMSCRIPTEN)
    # default to static libraries and no tests
    option(BUILD_SHARED_LIBS "Build shared library" OFF)
    option(BUILD_TESTING "Build tests" OFF)
else()
    # default to shared libraries and tests
    option(BUILD_SHARED_LIBS "Build shared library" ON)
    option(BUILD_TESTING "Build tests" ON)
endif()

# store whether we target shared lib
# so dependencies can be used statically if needed
set(IS_SHARED BUILD_SHARED_LIBS)

# Minimal standalone .wasm interface
cmake_dependent_option(VIPSTER_WEB "Build JS+WebGL GUI" ON "EMSCRIPTEN" OFF)
mark_as_advanced(VIPSTER_WEB)

# enable main GUI App if Qt is found
find_package(Qt6 CONFIG QUIET COMPONENTS Widgets Gui OpenGLWidgets)
if(NOT ${Qt6_FOUND})
    find_package(Qt5 5.10 CONFIG QUIET COMPONENTS Widgets Gui)
endif()

set(tmp ${VIPSTER_DESKTOP})
cmake_dependent_option(VIPSTER_DESKTOP "Build QT-based desktop app" ON "Qt6_FOUND OR Qt5_FOUND;NOT EMSCRIPTEN" OFF)
if(tmp AND NOT ${VIPSTER_DESKTOP})
    message(FATAL_ERROR "Requested to build the GUI, but no Qt installation was detected")
endif()

# embedded LAMMPS widget
cmake_dependent_option(VIPSTER_LAMMPS "Interactive LAMMPS in desktop app" OFF "VIPSTER_DESKTOP" OFF)

# Python widget and standalone library
find_package(Python3 3.6 QUIET COMPONENTS Interpreter Development)
cmake_dependent_option(VIPSTER_PYWIDGET "Python shell in desktop app" ${Python3_FOUND} "VIPSTER_DESKTOP" OFF)
option(VIPSTER_PYLIB "Standalone Python library" OFF)
mark_as_advanced(VIPSTER_PYLIB)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# include current dir so library can be accessed via "vipster/*.h"
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# add tinyexpr to the include path
include_directories("external/tinyexpr")

# Library, needed for everything else
file(GLOB_RECURSE PY_SOURCES "vipster/*.py.h" "vipster/*.py.cpp")
file(GLOB_RECURSE LIB_SOURCES "vipster/*.h" "vipster/*.tpp" "vipster/*.cpp" "external/tinyexpr.c")
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
        WINDOWS_EXPORT_ALL_SYMBOLS TRUE
    )
if(UNIX)
    # allow dynamic plugin loading
    target_link_libraries(libvipster PUBLIC dl)
    # on traditional unix-like platforms,
    # install files to appropriate locations
    if(BUILD_SHARED_LIBS AND NOT APPLE AND NOT VIPSTER_WEB)
        install(TARGETS libvipster
            EXPORT "Vipster"
            LIBRARY DESTINATION "lib"
            FRAMEWORK DESTINATION "lib"
            )
        install(DIRECTORY "vipster/"
            DESTINATION "include/vipster"
            FILES_MATCHING PATTERN "*.h"
            )
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h DESTINATION "include/vipster")
        target_include_directories(libvipster PUBLIC $<INSTALL_INTERFACE:include>)
        install(EXPORT "Vipster" DESTINATION "lib/cmake")
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
                "Please install ${NAME} or execute \"git submodule update --init external/${PATH}\".")
        else()
            message("Could not find suitable ${NAME} installation, using included version at external/${PATH}")
        endif()
    endif()
endfunction()

# find suitable json-library
find_package(nlohmann_json 3.10.5 CONFIG QUIET)
if(NOT nlohmann_json_FOUND)
    init_submodule("nlohmann/json" "nlohmann_json")
    set(JSON_BuildTests OFF CACHE INTERNAL "")
    set(JSON_Install OFF CACHE INTERNAL "")
    add_subdirectory(external/nlohmann_json EXCLUDE_FROM_ALL)
endif()
target_link_libraries(libvipster PRIVATE nlohmann_json::nlohmann_json)

# find suitable fmt-library
find_package(fmt 8.1.1 CONFIG QUIET)
if(NOT fmt_FOUND)
    init_submodule("{fmt}" "fmt")
    set(BUILD_SHARED_LIBS OFF)
    add_subdirectory(external/fmt EXCLUDE_FROM_ALL)
    set(BUILD_SHARED_LIBS IS_SHARED)
endif()
target_link_libraries(libvipster PRIVATE fmt::fmt)

# create python binding objects
if(VIPSTER_PYLIB OR VIPSTER_PYWIDGET)
    find_package(Python3 3.6 REQUIRED COMPONENTS Interpreter Development)
    find_package(pybind11 2.9.1 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        init_submodule("Pybind11" "pybind11")
        add_subdirectory(external/pybind11 EXCLUDE_FROM_ALL)
    endif()
    add_library(bindings OBJECT ${PY_SOURCES})
    target_link_libraries(bindings PUBLIC libvipster PRIVATE fmt::fmt pybind11::module)
    if(WIN32)
        add_definitions(-DHAVE_SNPRINTF)
    endif()
    if(NOT MSVC)
        target_compile_options(bindings PRIVATE "-fvisibility=hidden")
    endif()
endif()

if(VIPSTER_PYLIB)
    # standalone python bindings
    pybind11_add_module(pyvipster $<TARGET_OBJECTS:bindings> "vipster/main.py.cpp")
    set_target_properties(pyvipster PROPERTIES OUTPUT_NAME "vipster")
    target_link_libraries(pyvipster PRIVATE libvipster)
    if(VIPSTER_PYLIB OR (UNIX AND NOT APPLE))
        execute_process(COMMAND "${Python3_EXECUTABLE}" "-c"
"import sysconfig as s; print(s.get_path('platlib'))"
            OUTPUT_VARIABLE PYTHON_SITE
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        install(TARGETS pyvipster LIBRARY DESTINATION ${PYTHON_SITE})
        install(CODE  "execute_process(COMMAND \"${Python3_EXECUTABLE}\" \"${PROJECT_SOURCE_DIR}/setup.py\" \"egg_info\")")
        install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/vipster.egg-info" DESTINATION ${PYTHON_SITE})
    endif()
endif()

if(VIPSTER_DESKTOP)
    # Qt-based desktop application

    # load and setup qt
    find_package(Qt6 CONFIG COMPONENTS Widgets Gui OpenGLWidgets)
    if (NOT ${Qt6_FOUND})
        find_package(Qt5 CONFIG REQUIRED COMPONENTS Widgets Gui)
    endif()
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)

    # find suitable cli11 installation
    find_package(CLI11 2.1.2 CONFIG QUIET)
    if(NOT CLI11_FOUND)
        init_submodule("CLI11" "CLI11")
        set(CLI11_TESTING OFF CACHE INTERNAL "")
        add_subdirectory(external/CLI11 EXCLUDE_FROM_ALL)
    endif()

    # find OpenGL (needed if not in compiler's sysroot)
    set(OpenGL_GL_PREFERENCE GLVND)
    find_package(OpenGL)

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
    if(NOT MSVC)
        target_compile_options(vipster PRIVATE "-fvisibility=hidden")
    endif()

    # set rpath for local installations on linux
    if(UNIX AND BUILD_SHARED_LIBS AND NOT APPLE)
        set_target_properties(vipster PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
    endif()

    # fix name-clash on import library when using mingw/gcc
    if(WIN32)
        set_target_properties(vipster PROPERTIES IMPORT_SUFFIX "")
    endif()

    # link QT, CLI11 and vipster
    target_link_libraries(vipster PRIVATE
        Qt::Gui Qt::Widgets OpenGL::GL
        CLI11::CLI11
        fmt::fmt
        libvipster)
    if(${Qt6_FOUND})
        target_link_libraries(vipster PRIVATE Qt6::OpenGLWidgets)
    endif()
    target_include_directories(vipster PRIVATE "gui/common" "gui/qt")

    # python-widget:
    if(VIPSTER_PYWIDGET)
        add_definitions(-DUSE_PYTHON)
        target_sources(vipster PRIVATE ${QTPY_SOURCES} $<TARGET_OBJECTS:bindings>)
        target_link_libraries(vipster PRIVATE pybind11::embed Python3::Python)
    endif()

    # lammps-widget:
    if(VIPSTER_LAMMPS)
        set(USE_LAMMPS_BUILTIN ON)
        find_package(LAMMPS CONFIG QUIET)
        if(LAMMPS_FOUND)
            if(${LAMMPS_VERSION} GREATER_EQUAL 20201130)
                get_target_property(LMP_INC_DIR LAMMPS::lammps INTERFACE_INCLUDE_DIRECTORIES)
                if(EXISTS ${LMP_INC_DIR}/exceptions.h)
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
        install(FILES "util/com.github.sgsaenger.vipster.desktop" DESTINATION "share/applications")
        install(FILES "util/com.github.sgsaenger.vipster.appdata.xml" DESTINATION "share/metainfo")
    endif()
endif()

if(VIPSTER_WEB)
    # Emscripten-based web application
    set(EMCC_LINKER_FLAGS
        "--bind -s USE_WEBGL2=1 -s EXTRA_EXPORTED_RUNTIME_METHODS=['FS'] -s WASM=1 -s SINGLE_FILE=1 -s MODULARIZE=1 -s DISABLE_EXCEPTION_CATCHING=0 -s ALLOW_MEMORY_GROWTH=1"
        CACHE STRING
        "Configure how and what Emscripten bundles")
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

if(BUILD_TESTING)
    # enable testing via build system
    enable_testing()

    # Unit tests, so far only for library
    include(tests/unit_tests.cmake)

    # File I/O tests
    if(VIPSTER_DESKTOP AND NOT WIN32) # TODO: disable on Win due to path issues
        include(tests/file_tests.cmake)
    endif()
endif()
