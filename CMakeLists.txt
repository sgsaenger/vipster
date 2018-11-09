cmake_minimum_required(VERSION 3.9)

project(Vipster VERSION 1.14)

option(TESTS "Build tests" OFF)
option(DESKTOP "Build QT-based desktop app" ON)
include(CMakeDependentOption)
cmake_dependent_option(WEB "Build WASM-Module" OFF "NOT DESKTOP" OFF)
cmake_dependent_option(PYBIND "Python bindings" OFF "NOT WEB" OFF)

set(CMAKE_CXX_STANDARD 14)

include_directories(external)
include_directories(libvipster)
include_directories(python)
include_directories(vipster/qt)

# Library, needed for everything else
file(GLOB_RECURSE LIB_SOURCES "libvipster/*.cpp" "external/tinyexpr.c")
add_library(libvipster SHARED ${LIB_SOURCES})
set_target_properties(libvipster
    PROPERTIES
        OUTPUT_NAME "vipster"
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )

if(TESTS)
    # Tests, so far only for library
    file(GLOB TEST_SOURCES "tests/*.cpp")
    add_executable(test_lib ${TEST_SOURCES})
    target_link_libraries(test_lib libvipster)
endif()

if(PYBIND)
    # Python-bindings
    add_subdirectory(pybind11)
    file(GLOB PY_SOURCES "python/*.cpp")
    pybind11_add_module(pyvipster ${PY_SOURCES})
    set_target_properties(pyvipster PROPERTIES OUTPUT_NAME vipster)
    target_link_libraries(pyvipster PRIVATE libvipster)
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
     "vipster/qt/*.cpp"
     "vipster/qt/*.ui"
     "vipster/common/*.cpp"
     "vipster/resources/vipster.qrc"
     "vipster/resources/win.rc"
     )
add_executable(vipster ${QT_SOURCES})
target_link_libraries(vipster Qt5::Gui Qt5::Widgets libvipster)
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

file(GLOB WEB_SOURCES "vipster/common/*.cpp" "vipster/web/main.cpp")
add_executable(webvipster ${WEB_SOURCES})
target_link_libraries(webvipster libvipster)
endif()
