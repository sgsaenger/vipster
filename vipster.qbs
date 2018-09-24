import qbs 1.0

Project {
    name: "Vipster"
    minimumQbsVersion: "1.11.0"

    references: [
        "libvipster/libvipster.qbs",
        "vipster/QtVipster.qbs",
        "vipster/WebVipster.qbs",
        "python/python.qbs",
        "tests/tests.qbs",
        "dist/dist.qbs"
    ]

    // toggle emscripten-build
    property bool webBuild: false
    // create windows archive&&installer
    property bool winInstall: false
    // create osx install-.dmg
    property bool macInstall: false
    // enable python interface (TODO: not included in Win/OSX installers)
    property bool pythonBuild: false
    // specify non-default python executable (defaults to 'python' in $PATH)
    property string pythonName: "python"

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
