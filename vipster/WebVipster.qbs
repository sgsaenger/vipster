import qbs 1.0

CppApplication {
    name: "WebVipster"
    targetName: "vipster.js"
    condition: web
    Depends {name: "libvipster"}

    Rule {
        id: emscriptenHelper
        inputs: []
        multiplex: true
        outputFileTags: ["application"]
        outputArtifacts: [{
                filePath: product.destinationDirectory.concat('/').concat(
                              product.targetName.replace('js','wasm')),
                fileTags: ["application"]
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
                "resources/*vert",
                "default.json"]
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
        fileTagsFilter: ["application"]
        qbs.install: true
    }

    cpp.driverFlags: [
        "--embed-file", sourceDirectory.concat("/resources/atom.frag@atom.frag"),
        "--embed-file", sourceDirectory.concat("/resources/bond.frag@bond.frag"),
        "--embed-file", sourceDirectory.concat("/resources/cell.frag@cell.frag"),
        "--embed-file", sourceDirectory.concat("/resources/atom.vert@atom.vert"),
        "--embed-file", sourceDirectory.concat("/resources/bond.vert@bond.vert"),
        "--embed-file", sourceDirectory.concat("/resources/cell.vert@cell.vert"),
        "--embed-file", sourceDirectory.concat("/default.json@vipster.json")
    ]
}
