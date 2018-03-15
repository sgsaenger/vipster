import qbs 1.0

Product {
    name: "libvipster"
    targetName: "vipster"

    property bool isUnix: !qbs.targetOS.contains("windows") && !cpp.architecture.contains("wasm")

    type: "dynamiclibrary"
    Properties {
        condition: !isUnix
        type: "staticlibrary"
    }

    files: ["json.hpp",
            "libvipster.qmodel"]

    Group {
        name: "headers"
        files: ["*.h","*/*.h"]
        qbs.install: isUnix
        qbs.installDir: "include/vipster"
    }

    Group {
        name: "sources"
        files: ["*.cpp", "*/*.cpp"]
    }

    Group {
        name: "library"
        fileTagsFilter: product.type
        qbs.install: isUnix
        qbs.installDir: "lib"
    }

    // C++ settings (rest of project will inherit)
    Depends { name: "cpp"}
    cpp.cxxLanguageVersion: "c++14"
    Properties {
        condition: qbs.buildVariant=="profile"
        cpp.driverFlags: [
            "-fprofile-arcs",
            "-ftest-coverage"
        ]
        cpp.dynamicLibraries: ["gcov"]
    }
    Export {
        Depends { name: "cpp"}
        cpp.cxxLanguageVersion: "c++14"
        Properties {
            condition: qbs.buildVariant=="profile"
            cpp.driverFlags: [
                "-fprofile-arcs",
                "-ftest-coverage"
            ]
            cpp.dynamicLibraries: ["gcov"]
        }
        cpp.includePaths: [product.sourceDirectory]
    }
}
