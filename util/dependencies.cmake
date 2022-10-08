include(FetchContent)

set(FETCHCONTENT_QUIET OFF)

# library: tinyexpr
# purpose: general math expressions
# users:   library + GUI
# license: Zlib
FetchContent_Declare(
    tinyexpr
    GIT_REPOSITORY https://github.com/codeplea/tinyexpr
    GIT_TAG 373d2d3d5dc3344e3e2cb826717404064562a058
    GIT_SHALLOW TRUE
)

# library: nlohmann_json
# purpose: config file handling, file format support
# users:   library
# license: MIT
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json
    GIT_TAG v3.11.2
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS 3.11.2
)

# library: {fmt}
# purpose: string formatting
# users:   library + GUI
# license: MIT
# NOTE: LAMMPS has a dependency on fmt<9
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG 8.1.1
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS 8.1.1
)

# library: pybind11
# purpose: python bindings for library
# users:   pylib + GUI (optional)
# license: BSD
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.10.0
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS 2.10.0
)

# library: CLI11
# purpose: command-line parsing
# users:   GUI
# license: BSD
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11
    GIT_TAG v2.2.0
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS 2.2.0
)

# library: Catch2
# purpose: unit-testing
# users:   tests
# license: Boost
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2
    GIT_TAG v2.13.9
    GIT_SHALLOW TRUE
    FIND_PACKAGE_ARGS 2.13.9
)

# library: LAMMPS
# purpose: interactive molecular simulations
# users:   GUI
# license: GPL
FetchContent_Declare(
    lammps
    GIT_REPOSITORY https://github.com/lammps/lammps
    GIT_TAG stable_29Sep2021
    GIT_SHALLOW TRUE
)
