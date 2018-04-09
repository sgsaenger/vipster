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
            "qt/**"]

    // Install binary to install-root, under bin/ when it's a unix
    Group {
        name: "binary"
        fileTagsFilter: product.type
        qbs.install: !qbs.targetOS.contains("macos")
        Properties {
            condition: !qbs.targetOS.contains("windows")
            qbs.installDir: "bin"
        }
    }
    // Install OSX Bundle
    Group {
        fileTagsFilter: ["bundle.content"]
        qbs.install: true
        qbs.installDir: "."
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

    //TODO: OSX/Windows icon/packaging
//        Depends { name: "bundle" }
//        Depends { name: "ib" }
//        bundle.infoPlist: ({"CFBundleIconFile": "myapp"})
}
