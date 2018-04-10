import qbs 1.0

DynamicLibrary {
    name: "libvipster"
    targetName: "vipster"
    property bool isBundle: qbs.targetOS.contains("darwin") && bundle.isBundle

    files: ["json.hpp",
            "libvipster.qmodel"]

    Group {
        name: "headers"
        files: ["*.h","*/*.h"]
        qbs.install: !qbs.targetOS.contains("windows") && !qbs.targetOS.contains("macos") && !project.webBuild
        qbs.installDir: "include/vipster"
    }

    Group {
        name: "sources"
        files: ["*.cpp", "*/*.cpp"]
    }

    Group {
        name: "library"
        fileTagsFilter: isBundle ? ["bundle.content"] : product.type
        qbs.install: !project.webBuild
        qbs.installDir: isBundle ? "Vipster.app/Contents/Frameworks" : (qbs.targetOS.contains("windows") ? "" : "lib")
        qbs.installSourceBase: product.buildDirectory
    }

    Group {
        name: "distfiles"
        files: "default.json"
        qbs.install: !project.webBuild
        Properties {
            condition: !qbs.targetOS.contains("windows")
            qbs.installDir: {
                if (qbs.targetOS.contains("macos"))
                    return "Vipster.app/Contents/Resources"
                else
                    return "share/vipster"
            }
        }
    }

    // C++ settings (rest of project will inherit)
    Depends { name: "cpp"}
    cpp.cxxLanguageVersion: "c++14"
    //TODO: more intelligent wrapping. this is a rather bad idea!
    cpp.defines: [ "PREFIX="+project.prefix ]
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
        cpp.defines: [ "PREFIX="+project.prefix ]
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
