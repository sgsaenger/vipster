import qbs 1.0

Product {
    name: "WebVipster"
    targetName: "vipster.js"
    condition: project.webBuild
    type: "emapplication"
    Depends {name: "cpp"}
    Depends {name: "libvipster"}

    Rule {
        id: emscriptenHelper
        inputs: ["application"]
        multiplex: true
        outputFileTags: ["emapplication"]
        outputArtifacts: [{
                filePath: product.destinationDirectory.concat('/').concat(
                              product.targetName.replace('js','wasm')),
                fileTags: ["emapplication"]
            }]
        prepare: {
            var cmd = new Command('echo');
            cmd.silent = true;
            return [cmd];
        }
    }

    Group {
        name: "sources"
        files: ["common/*",
                "web/main.cpp"]
    }

    Group {
        // will be embedded in emscripten-product
        name: "resources"
        files: ["resources/*frag",
                "resources/*vert"]
    }

    Group {
        // HTML/JS-parts
        name: "page"
        prefix: "web/page/"
        files: ["index.html",
                "vipster_setup.js",
                "styles/styles.css"]
        qbs.install: true
    }

    Group {
        name: "binaries"
        fileTagsFilter: ["application", "emapplication"]
        qbs.install: true
    }

    cpp.driverFlags: [
        "--embed-file", sourceDirectory.concat("/resources/atom.frag@atom.frag"),
        "--embed-file", sourceDirectory.concat("/resources/bond.frag@bond.frag"),
        "--embed-file", sourceDirectory.concat("/resources/cell.frag@cell.frag"),
        "--embed-file", sourceDirectory.concat("/resources/select.frag@select.frag"),
        "--embed-file", sourceDirectory.concat("/resources/atom.vert@atom.vert"),
        "--embed-file", sourceDirectory.concat("/resources/bond.vert@bond.vert"),
        "--embed-file", sourceDirectory.concat("/resources/cell.vert@cell.vert"),
        "--embed-file", sourceDirectory.concat("/resources/select.vert@select.vert"),
        "--embed-file", sourceDirectory.concat("/resources/selection.frag@selection.frag"),
        "--embed-file", sourceDirectory.concat("/resources/selection.vert@selection.vert"),
        "--embed-file", sourceDirectory.concat("/../libvipster/default.json@vipster.json")
    ]
}
