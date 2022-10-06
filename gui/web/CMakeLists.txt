# Minimal standalone .wasm interface

file(GLOB WEB_SOURCES
    "../common/*.cpp"
    "main.cpp")

add_executable(webvipster ${WEB_SOURCES})

target_compile_definitions(webvipster PRIVATE WEBVIPSTER)

set_target_properties(webvipster PROPERTIES OUTPUT_NAME vipster)

target_include_directories(webvipster PRIVATE "../common")

target_link_libraries(webvipster PRIVATE libvipster)

file(GLOB WEB_RESOURCES "../resources/*.frag" "../resources/*.vert")
foreach(file ${WEB_RESOURCES})
    get_filename_component(filename ${file} NAME)
    set(EMCC_EMBED_FLAGS "${EMCC_EMBED_FLAGS}  --embed-file ${file}@${filename}")
endforeach(file)

target_link_options(webvipster PRIVATE
    "SHELL:--bind -s USE_WEBGL2=1 -s SINGLE_FILE=1 -s MODULARIZE=1 -s DISABLE_EXCEPTION_CATCHING=0 -s ALLOW_MEMORY_GROWTH=1"
    "SHELL:${EMCC_EMBED_FLAGS}"
)
