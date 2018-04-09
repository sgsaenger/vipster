import qbs 1.0

DynamicLibrary {
    name: "libvipster"
    targetName: "vipster"
    bundle.isBundle: false

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
            qbs.installDir: {
                if (qbs.targetOS.contains("macos"))
                    return "vipster.app/Contents/Frameworks"
                else
                    return  "lib"
            }
        }
    }

    Group {
        name: "distfiles"
        files: "default.json"
        qbs.install: !project.webBuild
        Properties {
            condition: !qbs.targetOS.contains("windows")
            qbs.installDir: {
                if (qbs.targetOS.contains("macos"))
                    return "vipster.app/Contents/Resources"
                else
                    return "share/vipster"
            }
        }
    }

    // C++ settings (rest of project will inherit)
    Depends { name: "cpp"}
    cpp.cxxLanguageVersion: "c++14"
    //TODO: more intelligent wrapping. this is a rather bad idea!
    cpp.defines: [ "PREFIX="+qbs.installRoot ]
    Properties {
        condition: qbs.targetOS.contains("macos")
        cpp.frameworks: ["CoreFoundation"]
        cpp.sonamePrefix: "@rpath"
    }
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
            condition: qbs.targetOS.contains("macos")
            cpp.frameworks: ["CoreFoundation"]
        }
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
