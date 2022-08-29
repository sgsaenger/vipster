file(GLOB_RECURSE TEST_SOURCES "tests/unit_tests/*.cpp")
add_executable(test_lib ${TEST_SOURCES})

# find suitable catch2-installation
find_package(Catch2 2.13.8 CONFIG QUIET)
if(NOT Catch2_FOUND)
    init_submodule("Catch2" "Catch2")
    set(CATCH_BUILD_TESTING OFF CACHE INTERNAL "")
    add_subdirectory(external/Catch2 EXCLUDE_FROM_ALL)
endif()
target_link_libraries(test_lib PRIVATE libvipster Catch2::Catch2)

add_test(NAME test_lib COMMAND $<TARGET_FILE:test_lib>)
