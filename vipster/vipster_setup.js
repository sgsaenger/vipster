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

function fillData() {
    var at = Module.getAtomIt(Module.curMol, Module.curStep);
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
    var cdm = document.getElementById("cellDim");
    cdm.value = Module.getCellDim(Module.curMol, Module.curStep, 0);
    var mat = Module.getCellVec(Module.curMol, Module.curStep);
    var vec = document.getElementById("cellVec");
    for(row=0; row<3; ++row){
        for(col=0; col<3; ++col){
            vec.rows[row+1].cells[col+1].innerText = mat[row][col];
        }
    }
}

function atomChanged(tgt) {
    var at = Module.getAtom(Module.curMol, Module.curStep,
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
    if (!isNaN(newVal)) {
        Module.setCellDim(Module.curMol, Module.curStep, newVal, 0);
    }
    Module.updateView();
    fillData();
}

function cellEnabled(val) {
    Module.enableCell(Module.curMol, Module.curStep, val);
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
    Module.setCellVec(Module.curMol, Module.curStep, vec);
    Module.updateView();
    fillData();
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
    document.getElementById('stepCur').innerHTML = 1;
    slider.value = 0;
    slider.max = nstep - 1;
    Module.curStep = 0;
    Module.setStep(i, 0);
    fillData();
}

function setStep(e) {
    var i = parseInt(e.target.value);
    document.getElementById('stepCur').innerHTML = i + 1;
    Module.curStep = i;
    Module.setStep(Module.curMol, i);
    fillData();
}
