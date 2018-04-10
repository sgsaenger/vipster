import qbs 1.0

QtGuiApplication {
    name: "QtVipster"
    targetName: "vipster"
    condition: !project.webBuild
    Depends { name: "libvipster" }
    Depends { name: "Qt"; submodules: ["core", "widgets", "opengl"] }
    cpp.includePaths: [product.sourceDirectory]

    property bool isBundle: qbs.targetOS.contains("darwin") && bundle.isBundle

    files: ["common/*",
            "resources/vipster.qrc",
            "qt/**"]

    Group {
        fileTagsFilter: isBundle ? ["bundle.content"] : ["application"]
        qbs.install: true
        qbs.installDir: isBundle ? "." : (qbs.targetOS.contains("windows") ? "" : "bin")
        qbs.installSourceBase: product.buildDirectory
    }

    Group {
        name: "WinIcon"
        files: ["resources/vipster.ico",
                "resources/win.rc"]
    }
    Group {
        name: "OSXIcon"
        files: ["resources/vipster.icns"]
    }

    Properties{
        condition: qbs.targetOS.contains("macos")
        bundle.infoPlist: ({"CFBundleIconFile": "vipster"})
        cpp.useRPaths: true
        cpp.rpaths: ["@loader_path/../Frameworks"]
    }
}
