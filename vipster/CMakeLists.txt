# libvipster
# main library providing containers and file format support

# Dependencies
FetchContent_MakeAvailable(nlohmann_json fmt tinyexpr)

file(GLOB_RECURSE LIB_HEADERS "*.h")
configure_file("version.h.in" "version.h")

file(GLOB_RECURSE LIB_SOURCES "*.tpp" "*.cpp" "${tinyexpr_SOURCE_DIR}/tinyexpr.c")

add_library(libvipster ${LIB_HEADERS} ${LIB_SOURCES})

set_target_properties(libvipster
    PROPERTIES
        FRAMEWORK TRUE
        OUTPUT_NAME "vipster"
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        WINDOWS_EXPORT_ALL_SYMBOLS TRUE
)

target_include_directories(libvipster PUBLIC
    $<BUILD_INTERFACE:${tinyexpr_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(libvipster PRIVATE
    ${CMAKE_DL_LIBS}
    nlohmann_json::nlohmann_json
    fmt::fmt-header-only
)

if(NOT EMSCRIPTEN)
    install(TARGETS libvipster
        EXPORT "Vipster"
        LIBRARY DESTINATION "lib"
        FRAMEWORK DESTINATION "lib"
    )
    install(DIRECTORY "."
        DESTINATION "include/vipster"
        FILES_MATCHING PATTERN "*.h"
    )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.h DESTINATION "include/vipster")
    install(EXPORT "Vipster" DESTINATION "lib/cmake")
endif()
