import qbs 1.0

Product {
    name: "libvipster"
    targetName: "vipster"

    type: "dynamiclibrary"

    files: ["json.hpp",
            "libvipster.qmodel"]

    Group {
        name: "headers"
        files: ["*.h","*/*.h"]
        qbs.install: !qbs.targetOS.contains("windows") && !project.webBuild
        qbs.installDir: "include/vipster"
    }

    Group {
        name: "sources"
        files: ["*.cpp", "*/*.cpp"]
    }

    Group {
        name: "library"
        fileTagsFilter: product.type
        qbs.install: !project.webBuild
        Properties {
            condition: !qbs.targetOS.contains("windows")
            qbs.installDir: "lib"
        }
    }

    Group {
        name: "distfiles"
        files: "default.json"
        qbs.install: !project.webBuild
        Properties {
            condition: !qbs.targetOS.contains("windows")
            qbs.installDir: "share/vipster"
        }
    }

    // C++ settings (rest of project will inherit)
    Depends { name: "cpp"}
    cpp.cxxLanguageVersion: "c++14"
    cpp.defines: [ "PREFIX="+qbs.installRoot ]
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
        cpp.defines: [ "PREFIX="+qbs.installRoot ]
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
