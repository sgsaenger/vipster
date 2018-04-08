import qbs 1.0

Project {
    name: "Vipster"

    references: [
        "libvipster/libvipster.qbs",
        "vipster/QtVipster.qbs",
        "vipster/WebVipster.qbs",
        "python/python.qbs",
        "tests/tests.qbs",
        "dist/dist.qbs"
    ]

    property bool webBuild: false
    property bool winInstall: false
    property bool pythonBuild: false
    property string pythonName: ""

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
