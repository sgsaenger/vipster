import qbs 1.0

QtGuiApplication {
    name: "QtVipster"
    targetName: qbs.targetOS.contains("macos") ? "Vipster" : "vipster"
    condition: !project.webBuild
    consoleApplication: false
    Depends {
        name: "libvipster"
    }
    Depends {
        name: "Qt"
        submodules: ["core", "widgets", "opengl"]
    }
    cpp.includePaths: [product.sourceDirectory]

    property bool isBundle: qbs.targetOS.contains("darwin") && bundle.isBundle

    files: ["common/*", "resources/vipster.qrc", "qt/**"]

    Group {
        fileTagsFilter: isBundle ? ["bundle.content"] : ["application"]
        qbs.install: true
        qbs.installDir: isBundle ? "." : (qbs.targetOS.contains(
                                              "windows") ? "" : "bin")
        qbs.installSourceBase: product.buildDirectory
    }

    Group {
        name: "WinIcon"
        files: ["resources/vipster.ico", "resources/win.rc"]
    }
    Group {
        name: "OSXIcon"
        files: ["resources/vipster.icns"]
    }

    Properties {
        condition: qbs.targetOS.contains("macos")
        bundle.infoPlist: ({
                               CFBundleIconFile: "vipster"
                           })
        cpp.useRPaths: true
        cpp.rpaths: ["@loader_path/../Frameworks"]
    }

    Group {
        name: "LinuxIcon"
        files: ["resources/vipster.png"]
        qbs.install: !qbs.targetOS.contains("windows")
                     && !qbs.targetOS.contains("macos")
        qbs.installDir: "share/icons/hicolor/128x128/apps"
    }

    Group {
        name: "DesktopFile"
        files: ["resources/vipster.desktop"]
        qbs.install: !qbs.targetOS.contains("windows")
                     && !qbs.targetOS.contains("macos")
        qbs.installDir: "share/applications"
    }

    Group {
        name: "AppstreamFile"
        files: ["resources/vipster.appdata.xml"]
        qbs.install: !qbs.targetOS.contains("windows")
                     && !qbs.targetOS.contains("macos")
        qbs.installDir: "share/metainfo"
    }
}
