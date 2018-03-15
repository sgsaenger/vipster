import qbs 1.0

CppApplication {
    name: "test_libvipster"
    condition: !cpp.architecture.contains("wasm")
    files: ["catch.hpp",
            "*.cpp"]
    Properties {
        condition: qbs.targetOS.contains("windows")
        cpp.defines: "DO_NOT_USE_WMAIN"
        cpp.driverFlags: [
            "-static",
            "-static-libstdc++",
            "-static-libgcc"
        ]
    }
    Depends { name: "libvipster" }
}
