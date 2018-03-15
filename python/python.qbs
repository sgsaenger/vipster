import qbs 1.0

DynamicLibrary {
    name: "PyVipster"
    builtByDefault: false
    files: "pyvipster.cpp"
    Depends { name: "libvipster" }
}
