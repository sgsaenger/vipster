import qbs 1.0

QtGuiApplication {
    name: "QtVipster"
    targetName: "vipster"
    condition: !web
    Depends { name: "libvipster" }
    Depends { name: "Qt"; submodules: ["core", "widgets", "opengl"] }
    cpp.includePaths: [product.sourceDirectory]

    files: ["common/*",
            "resources/vipster.qrc",
            "resources/vipster.icns",
            "qt/*"]

    Group {
        name: "binary"
        fileTagsFilter: product.type
        qbs.install: true
        qbs.installDir: "bin"
    }

    Group {
        name: "distfiles"
        files: "default.json"
        qbs.install: true
        Properties {
            condition: qbs.targetOS.contains("windows")
            qbs.installDir: "share/vipster"
        }
    }

    //TODO: OSX/Windows icon/packaging
//        Depends { name: "bundle" }
//        Depends { name: "wix" }
//        Depends { name: "ib" }
//        bundle.infoPlist: ({"CFBundleIconFile": "myapp"})
}
