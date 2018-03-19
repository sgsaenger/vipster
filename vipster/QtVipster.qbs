import qbs 1.0

QtGuiApplication {
    name: "QtVipster"
    targetName: "vipster"
    condition: !project.webBuild
    Depends { name: "libvipster" }
    Depends { name: "Qt"; submodules: ["core", "widgets", "opengl"] }
    cpp.includePaths: [product.sourceDirectory]

    files: ["common/*",
            "resources/vipster.qrc",
            "qt/*"]

    Group {
        name: "binary"
        fileTagsFilter: product.type
        qbs.install: true
        Properties {
            condition: !qbs.targetOS.contains("windows")
            qbs.installDir: "bin"
        }
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


    //TODO: OSX/Windows icon/packaging
//        Depends { name: "bundle" }
//        Depends { name: "ib" }
//        bundle.infoPlist: ({"CFBundleIconFile": "myapp"})
}
