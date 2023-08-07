file(GLOB_RECURSE TEST_SOURCES "unit_tests/*.cpp")
add_executable(test_lib ${TEST_SOURCES})

FetchContent_MakeAvailable(Catch2)
target_link_libraries(test_lib PRIVATE libvipster Catch2::Catch2)

add_test(NAME test_lib COMMAND $<TARGET_FILE:test_lib>)
set_property(TEST test_lib PROPERTY ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${CMAKE_BINARY_DIR}/vipster")
