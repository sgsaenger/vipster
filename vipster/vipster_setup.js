var Change = {
    atoms: 1,
    cell: 2,
    fmt: 4,
    //
    kpoints: 32,
};
Change.step =  Change.atoms | Change.cell | Change.fmt;
Change.mol = Change.kpoints;

var Module = {
    preRun: [],
    postRun: [],
    curMol: 0,
    curStep: 0,
    canvas: (function () {
        var canvas = document.getElementById('canvas');
        canvas.addEventListener("webglcontextlost", function (e) {
            alert('WebGL context lost. You will need to reload the page.');
            e.preventDefault();
        }, false);
        canvas.width = canvas.clientWidth;
        canvas.height = canvas.clientHeight;
        return canvas;
    })()
}

function update(change) {
    if (change & (Change.atoms | Change.cell))
        fillAtoms();
    if (change & Change.cell)
        fillCell();
}

function fillAtoms() {
    var fmt = parseInt(document.getElementById('atFmtSel').value);
    var at = Module.getAtomIt(Module.curMol, Module.curStep, fmt);
    var nat = Module.getNAtoms(Module.curMol, Module.curStep);
    var atList = document.getElementById('atList');
    atList.innerHTML = "<tr><th>Type</th><th>X</th><th>Y</th><th>Z</th></tr>";
    for (i = 0; i < nat; ++i) {
        atList.innerHTML += "<tr data-idx=" + i + "> \
<td contenteditable data-idx='name'>" + at.name + "</td> \
<td contenteditable data-idx='0'>" + at.coord[0] + "</td> \
<td contenteditable data-idx='1'>" + at.coord[1] + "</td> \
<td contenteditable data-idx='2'>" + at.coord[2] + "</td></tr>";
        at.increment();
    }
}

function fillCell() {
    var cdm = document.getElementById("cellDim");
    var fmt = parseInt(document.getElementById("cdmFmtSel").value);
    cdm.value = Module.getCellDim(Module.curMol, Module.curStep, fmt);
    var mat = Module.getCellVec(Module.curMol, Module.curStep);
    var vec = document.getElementById("cellVec");
    for(row=0; row<3; ++row){
        for(col=0; col<3; ++col){
            vec.rows[row+1].cells[col+1].innerText = mat[row][col];
        }
    }
}

function atomChanged(tgt) {
    var fmt = parseInt(document.getElementById('atFmtSel').value);
    var at = Module.getAtom(Module.curMol, Module.curStep, fmt,
                            parseInt(tgt.parentElement.dataset.idx));
    if (tgt.dataset.idx === 'name') {
        at.name = tgt.innerText;
    } else {
        var dir = parseInt(tgt.dataset.idx);
        var newVal = parseFloat(tgt.innerText);
        if (!isNaN(newVal)) {
            var coord = at.coord;
            coord[dir] = newVal;
            at.coord = coord;
        }
    }
    Module.updateView();
}

function cellDimChanged(tgt) {
    var newVal = parseFloat(tgt.value);
    if (isNaN(newVal))
        return;
    var fmt = parseInt(document.getElementById("cdmFmtSel").value);
    var scale = document.getElementById("cellScale").checked;
    var trajec = document.getElementById("cellToMol").checked;
    if (trajec) {
        for (var i=0;i<Module.getMolNStep(Module.curMol);++i) {
            Module.setCellDim(Module.curMol, i, newVal, fmt, scale);
        }
    } else {
        Module.setCellDim(Module.curMol, Module.curStep, newVal, fmt, scale);
    }
    Module.updateView();
    update(Change.cell);
}

function cellEnabled(val) {
    var trajec = document.getElementById("cellToMol").checked;
    if (trajec) {
        for (var i=0;i<Module.getMolNStep(Module.curMol);++i) {
            Module.enableCell(Module.curMol, i, val);
        }
    } else {
        Module.enableCell(Module.curMol, Module.curStep, val);
    }
    Module.updateView();
}

function cellVecChanged(tgt) {
    var newVal = parseFloat(tgt.innerText);
    if (isNaN(newVal))
        return;
    var col = tgt.dataset.idx;
    var row = tgt.parentElement.dataset.idx;
    var vec = Module.getCellVec(Module.curMol, Module.curStep);
    vec[row][col] = newVal;
    var trajec = document.getElementById("cellToMol").checked;
    var scale = document.getElementById("cellScale").checked;
    if (trajec) {
        for (var i=0;i<Module.getMolNStep(Module.curMol);++i) {
            Module.setCellVec(Module.curMol, i, vec, scale);
        }
    } else {
        Module.setCellVec(Module.curMol, Module.curStep, vec, scale);
    }
    Module.updateView();
    update(Change.cell);
}

function readFile(e) {
    var file = document.getElementById('upfile').files[0];
    var molList = document.getElementById('molList');
    var reader = new FileReader();
    reader.onload = function (e) {
        Module.FS_createDataFile("/tmp", "test.file", e.target.result, true);
        Module.readFile("/tmp/test.file", file.name,
                        parseInt(document.getElementById('uptype').value));
        Module.FS_unlink("/tmp/test.file");
        molList.innerHTML += '<li onclick="setMol(getMolListIdx(event.target))">'
                + file.name + "</li>";
        setMol(molList.childElementCount - 1);
    }
    reader.readAsText(file);
}

function setMult(e) {
    var x = document.getElementById('xmult').value;
    var y = document.getElementById('ymult').value;
    var z = document.getElementById('zmult').value;
    Module.setMult(parseInt(x), parseInt(y), parseInt(z));
}

function getMolListIdx(li) {
    var molList = li.parentElement.children;
    for (i = 0; i < li.parentElement.childElementCount; i++) {
        if (li === molList[i])
            break;
    }
    return i;
}

function setMol(i) {
    Module.curMol = i;
    var nstep = Module.getMolNStep(i);
    var slider = document.getElementById('stepSlider');
    document.getElementById('stepMax').innerHTML = nstep;
    slider.max = nstep - 1;
    slider.value = nstep - 1;
    setStep(nstep - 1);
}

function setStep(i) {
    document.getElementById('stepCur').innerHTML = i + 1;
    Module.curStep = i;
    Module.setStep(Module.curMol, i);
    document.getElementById('atFmtSel').value = Module.getFmt(Module.curMol, Module.curStep);
    update(Change.step);
}

function addParser(idx, name) {
    document.getElementById('uptype').innerHTML +=
            '<option value='+idx+'>'+UTF8ToString(name)+'</option>';
}
