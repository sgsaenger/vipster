# Full-featured Qt-based GUI

# setup qt
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# get source files and split off optional sources
file(GLOB_RECURSE QT_SOURCES
     "*.h"
     "*.cpp"
     "*.ui"
     "../common/*.h"
     "../common/*.cpp"
     "../resources/vipster.qrc"
     )
file(GLOB_RECURSE QTPY_SOURCES "*.py.*")
list(REMOVE_ITEM QT_SOURCES ${QTPY_SOURCES})
file(GLOB_RECURSE QTLMP_SOURCES "*.lmp.*")
list(REMOVE_ITEM QT_SOURCES ${QTLMP_SOURCES})

# icon-files:
if(APPLE)
    set(MACOSX_BUNDLE_ICON_FILE vipster.icns)
    set(ICON_SOURCES "${PROJECT_SOURCE_DIR}/util/vipster.icns")
    set_source_files_properties(${ICON_SOURCES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
elseif(WIN32)
    set(ICON_SOURCES "${PROJECT_SOURCE_DIR}/util/win.rc")
endif()

# executable
qt_add_executable(vipster WIN32 MACOSX_BUNDLE ${QT_SOURCES} ${ICON_SOURCES})
if(NOT MSVC)
    target_compile_options(vipster PRIVATE "-fvisibility=hidden")
endif()

# on traditional unix-like platforms,
# install files to appropriate locations
if(UNIX AND NOT APPLE)
    install(TARGETS vipster RUNTIME DESTINATION "bin")
    install(FILES "${PROJECT_SOURCE_DIR}/util/vipster.png" DESTINATION "share/icons/hicolor/128x128/apps")
    install(FILES "${PROJECT_SOURCE_DIR}/util/com.github.sgsaenger.vipster.desktop" DESTINATION "share/applications")
    install(FILES "${PROJECT_SOURCE_DIR}/util/com.github.sgsaenger.vipster.appdata.xml" DESTINATION "share/metainfo")
endif()
# set rpath for local installations on linux
if(UNIX AND BUILD_SHARED_LIBS AND NOT APPLE)
    set_target_properties(vipster PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
endif()

# fix name-clash on import library when using mingw/gcc
if(MINGW)
    set_target_properties(vipster PROPERTIES IMPORT_SUFFIX "")
endif()

# dependencies
FetchContent_MakeAvailable(CLI11 fmt)
target_link_libraries(vipster PRIVATE
    Qt::Gui Qt::Widgets
    CLI11::CLI11
    fmt::fmt-header-only
    libvipster)
if(${Qt6_FOUND})
    target_link_libraries(vipster PRIVATE Qt6::OpenGLWidgets)
endif()

# ensure availability of OpenGL
if(EMSCRIPTEN)
    target_link_options(vipster PRIVATE "SHELL:-s FULL_ES3=1 -s USE_WEBGL2=1")
else()
    set(OpenGL_GL_PREFERENCE GLVND)
    find_package(OpenGL)
    target_link_libraries(vipster PRIVATE OpenGL::GL)
endif()

target_include_directories(vipster PRIVATE "../common" ".")

# python-widget:
if(VIPSTER_PYWIDGET)
    FetchContent_MakeAvailable(pybind11)
    add_definitions(-DUSE_PYTHON)
    target_sources(vipster PRIVATE ${QTPY_SOURCES})
    target_link_libraries(vipster PRIVATE bindings pybind11::embed Python3::Python)
endif()

# lammps-widget:
if(VIPSTER_LAMMPS)
    set(USE_LAMMPS_BUILTIN ON)
    #find_package(LAMMPS CONFIG QUIET)
    #if(LAMMPS_FOUND)
    #    if(${LAMMPS_VERSION} GREATER_EQUAL 20201130)
    #        get_target_property(LMP_INC_DIR LAMMPS::lammps INTERFACE_INCLUDE_DIRECTORIES)
    #        if(EXISTS ${LMP_INC_DIR}/exceptions.h)
    #            message("Found LAMMPS (suitable version \"${LAMMPS_VERSION}\" at ${LMP_INC_DIR})")
    #            set(USE_LAMMPS_BUILTIN OFF)
    #            # link against imported target
    #            target_link_libraries(vipster PRIVATE LAMMPS::lammps)
    #        else()
    #            message("Found incompatible LAMMPS (suitable version \"${LAMMPS_VERSION}\" at ${LMP_INC_DIR} has exceptions disabled)")
    #        endif()
    #    else()
    #        message("Found incompatible LAMMPS (version  \"${LAMMPS_VERSION}\" too low)")
    #    endif()
    #endif()
    if(USE_LAMMPS_BUILTIN)
        FetchContent_GetProperties(lammps)
        if(NOT LAMMPS_POPULATED)
            FetchContent_Populate(lammps)
        endif()
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
        message(STATUS "Lammps: ${lammps_SOURCE_DIR} ${lammps_BINARY_DIR}")
        add_subdirectory(${lammps_SOURCE_DIR}/cmake ${lammps_BINARY_DIR} EXCLUDE_FROM_ALL)
        set(BUILD_SHARED_LIBS IS_SHARED)
        # link against own built target (LAMMPS::lammps may be mis-configured by find_package)
        target_link_libraries(vipster PRIVATE lammps)
    endif()
    # enable widget in vipster
    add_definitions(-DUSE_LAMMPS)
    target_sources(vipster PRIVATE ${QTLMP_SOURCES})
endif()
