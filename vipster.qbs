import qbs 1.0

Project {
    name: "Vipster"

    references: [
        "libvipster/libvipster.qbs",
//        "python/python.qbs",
        "tests/tests.qbs"
    ]

    SubProject {
        filePath: "vipster/frontend.qbs"
        inheritProperties: true
        Properties {
            web: profile.contains("mscript")
        }
    }

    Profile {
        name: "emscripten"
        cpp.compilerName: "emcc"
        cpp.cxxCompilerName: "em++"
        cpp.architecture: "wasm"
        cpp.driverFlags: ["--bind",
            "-s", "ERROR_ON_UNDEFINED_SYMBOLS=1",
            "-s", "USE_WEBGL2=1",
            "-s", "WASM=1",
            "-s", "DISABLE_EXCEPTION_CATCHING=0",
            "-O3"
        ]
        cpp.compilerIncludePaths: undefined
        cpp.compilerFrameworkPaths: undefined
        cpp.compilerLibraryPaths: undefined
    }
}
