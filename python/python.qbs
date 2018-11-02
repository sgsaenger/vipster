import qbs 1.0
import qbs.Probes as Probes
import qbs.Process

DynamicLibrary {
    name: "PyVipster"
    targetName: "vipster" + pySettings.libSuffix
    condition: hasPython.found && !project.webBuild
    builtByDefault: project.pythonBuild
    files: ["*.h", "*.cpp", "*.pypp"]
    Depends { name: "libvipster" }

    Group{
        name: "target"
        fileTagsFilter: product.type
        qbs.install: !qbs.targetOS.contains("windows") && !qbs.targetOS.contains("macos")
        qbs.installDir: pySettings.instPath
    }

    // fix unix library-naming convention:
    cpp.staticLibraryPrefix: ''
    cpp.dynamicLibraryPrefix: ''
    cpp.loadableModulePrefix: ''

    Probes.BinaryProbe {
        id: hasPython
        names: project.pythonName
    }

    Probe{
        id: pySettings
        condition: hasPython.found
        property string incPath
        property string instPath
        property string libSuffix
        configure: {
            var proc = new Process();
            proc.exec(hasPython.filePath, ["-c",
                      "import sysconfig as sc;"+
                      "print(sc.get_path('include'));"+
                      "print(sc.get_path('platlib', vars={'base':'','platbase':'','userbase':''}));"+
                      "print(sc.get_config_var('SOABI'));"
                      ]);
            incPath = proc.readLine();
            instPath = proc.readLine();
            instPath = instPath.substr(1);
            libSuffix = proc.readLine();
            if(libSuffix=="None"){
                libSuffix = undefined;
            }else{
                libSuffix = "."+libSuffix;
            }
        }
    }

    cpp.includePaths: [pySettings.incPath]
}
