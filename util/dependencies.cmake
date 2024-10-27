include(FetchContent)

# library: tinyexpr
# purpose: general math expressions
# users:   library + GUI
# license: Zlib
set(tinyexpr_TAG 373d2d3d5dc3344e3e2cb826717404064562a058)
set(tinyexpr_REPOSITORY https://github.com/codeplea/tinyexpr)

# library: nlohmann_json
# purpose: config file handling, file format support
# users:   library
# license: MIT
set(nlohmann_json_VERSION 3.11.2)
set(nlohmann_json_TAG v3.11.2)
set(nlohmann_json_REPOSITORY https://github.com/nlohmann/json)

# library: {fmt}
# purpose: string formatting
# users:   library + GUI
# license: MIT
set(fmt_VERSION 10.1.1)
set(fmt_TAG 10.1.1)
set(fmt_REPOSITORY https://github.com/fmtlib/fmt)

# library: pybind11
# purpose: python bindings for library
# users:   pylib + GUI (optional)
# license: BSD
set(pybind11_VERSION 2.10.0)
set(pybind11_TAG v2.10.0)
set(pybind11_REPOSITORY https://github.com/pybind/pybind11)

# library: CLI11
# purpose: command-line parsing
# users:   GUI
# license: BSD
set(CLI11_VERSION 2.2.0)
set(CLI11_TAG v2.2.0)
set(CLI11_REPOSITORY https://github.com/CLIUtils/CLI11)

# library: Catch2
# purpose: unit-testing
# users:   tests
# license: Boost
set(Catch2_VERSION 2.13.9)
set(Catch2_TAG v2.13.9)
set(Catch2_REPOSITORY https://github.com/catchorg/Catch2)

# library: LAMMPS
# purpose: interactive molecular simulations
# users:   GUI
# license: GPL
set(lammps_TAG stable_2Aug2023)
set(lammps_REPOSITORY https://github.com/lammps/lammps)

set(FETCHCONTENT_QUIET OFF)

foreach (dep IN ITEMS nlohmann_json fmt pybind11 CLI11 Catch2 lammps)
    message(STATUS "Registering ${dep} from ${${dep}_REPOSITORY} with Tag ${${dep}_TAG}")
    FetchContent_Declare(
        ${dep}
        GIT_REPOSITORY ${${dep}_REPOSITORY}
        GIT_TAG ${${dep}_TAG}
        GIT_SHALLOW TRUE
        FIND_PACKAGE_ARGS ${${dep}_VERSION}
    )
endforeach()

if(VIPSTER_DOWNLOAD_DEPENDENCIES)
    function(find_or_download_package PACKAGE)
        FetchContent_MakeAvailable(${PACKAGE})
    endfunction()

    # special dependency: tinyexpr cannot use GIT_SHALLOW
    FetchContent_Declare(
        tinyexpr
        GIT_REPOSITORY ${tinyexpr_REPOSITORY}
        GIT_TAG ${tinyexpr_TAG}
    )
    FetchContent_MakeAvailable(tinyexpr)
else()
    function(find_or_download_package PACKAGE)
        find_package(${PACKAGE} ${${PACKAGE}_VERSION} REQUIRED)
    endfunction()

    # special dependency: tinyexpr cannot be consumed find_package, provide path manually
    if (NOT DEFINED tinyexpr_SOURCE_DIR)
        message(FATAL_ERROR "If VIPSTER_DOWNLOAD_DEPENDENCIES is disabled, tinyexpr_SOURCE_DIR must be provided")
    endif()
endif()
