# Python bindings for libvipster

# actual bindings
# provided as an object library for further cmake consumption
file(GLOB_RECURSE PY_SOURCES "*.py.h" "*.py.cpp")
list(REMOVE_ITEM PY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.py.cpp")
add_library(bindings OBJECT ${PY_SOURCES})

set_target_properties(bindings PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE # make linkable into .so
    EXCLUDE_FROM_ALL TRUE # only build when required by another target
)

# dependencies
FetchContent_MakeAvailable(pybind11 fmt)
target_link_libraries(bindings PUBLIC libvipster pybind11::pybind11 PRIVATE fmt::fmt-header-only)

target_include_directories(bindings PUBLIC ".")

if(WIN32)
    # maybe no longer needed since Py 3.9
    add_definitions(-DHAVE_SNPRINTF)
endif()
if(NOT MSVC)
    # reduce library size
    target_compile_options(bindings PRIVATE "-fvisibility=hidden")
endif()

# standalone python library target
if(VIPSTER_PYLIB)
    pybind11_add_module(pyvipster "main.py.cpp")

    set_target_properties(pyvipster PROPERTIES OUTPUT_NAME "vipster")

    target_link_libraries(pyvipster PRIVATE bindings pybind11::module)

    execute_process(COMMAND "${Python3_EXECUTABLE}" "-c"
"import sysconfig as s; print(s.get_path('platlib'))"
        OUTPUT_VARIABLE PYTHON_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    install(TARGETS pyvipster LIBRARY DESTINATION ${PYTHON_SITE})
    install(CODE  "execute_process(COMMAND \"${Python3_EXECUTABLE}\" \"${PROJECT_SOURCE_DIR}/setup.py\" \"egg_info\")")
    install(DIRECTORY "${PROJECT_BINARY_DIR}/vipster.egg-info" DESTINATION ${PYTHON_SITE})
endif()
