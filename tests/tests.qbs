import qbs 1.0

CppApplication {
    name: "test_libvipster"
    condition: !project.webBuild
    files: ["catch.hpp",
            "*.cpp"]
    Properties {
        condition: qbs.targetOS.contains("windows")
        cpp.defines: "DO_NOT_USE_WMAIN"
    }
    Depends { name: "libvipster" }
    // install binary on windows because qbs does not add vipster.dll to search path
    property bool qbsfix: false
    Group {
        name: binary
        fileTagsFilter: product.type
        qbs.install: qbsfix
    }
}
