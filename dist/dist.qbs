import qbs 1.0
import qbs.FileInfo

Project {
    property string baseName:{
            var base = "Vipster-";
            if (qbs.targetOS.contains("windows")){
                base += "Win-";
            } else if (qbs.targetOS.contains("macos")){
                base += "OSX-";
            }
            if (qbs.debugInformation){
                base += "debug-";
            }
            if (qbs.architecture === "x86_64") {
                base += "64"
            } else {
                base += "32"
            }
            return base;
        }
    Product {
        name: "Qt DLLs"
        condition: project.winInstall
        Depends {name: "Qt.core"}
        property string suffix: {
            if (qbs.debugInformation) {
                return "d.dll";
            } else {
                return ".dll";
            }
        }
        Group {
            prefix: Qt.core.binPath + '/'
            files: [
                    "Qt5Core" + suffix,
                    "Qt5Gui" + suffix,
                    "Qt5Widgets" + suffix,
                    "libgcc_s_dw2-1.dll",
                    "libwinpthread-1.dll",
                    "libstdc++-6.dll",
            ]
            qbs.install: true
        }
        Group {
            prefix: Qt.core.pluginPath + "/platforms/"
            files: ["qwindows"+suffix]
            qbs.install: true
            qbs.installDir: "platforms"
        }
    }

    InstallPackage {
        name: "winArchive"
        condition: project.winInstall
        builtByDefault: true
        targetName: project.baseName
        archiver.type: "zip"
        Depends {name: "libvipster"}
        Depends {name: "QtVipster"}
        Depends {name: "Qt DLLs"}
        destinationDirectory: project.buildDirectory
    }

    NSISSetup {
        name: "winSetup"
        condition: project.winInstall
        targetName: project.baseName + "-install"
        Depends {name: "libvipster"}
        Depends {name: "QtVipster"}
        Depends {name: "Qt DLLs"}
        files: ["win.nsi"]
        nsis.defines: [
            "buildDirectory=" + FileInfo.toWindowsSeparators(qbs.installRoot),
        ]
        nsis.compressor: "lzma"
        destinationDirectory: project.buildDirectory
    }

    AppleApplicationDiskImage {
        name: "macSetup"
        condition: project.macInstall
        targetName: "Vipster"
        Depends {name: "QtVipster"}
        Depends {name: "libvipster"}
        version: project.version
        files: ["../vipster/resources/vipster.icns"]
        sourceBase: ""
        dmg.volumeName: "Vipster"
        dmg.iconSize: 128
        dmg.windowWidth: 640
        dmg.windowHeight: 502
        dmg.iconPositions: [
            {"path": "Applications", "x": 128, "y": 128},
            {"path": "Vipster.app", "x": 256, "y": 128},

        ]
        Group {
         fileTagsFilter: product.type
         qbs.install:true
        }
    }
}
