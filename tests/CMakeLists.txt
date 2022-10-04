# Unit tests, so far only for library
include(unit_tests.cmake)

# File I/O tests
if(VIPSTER_DESKTOP AND NOT WIN32) # TODO: disabled on Win due to path issues
    include(file_tests.cmake)
endif()
