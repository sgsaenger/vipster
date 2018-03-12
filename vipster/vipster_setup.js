const VERBOSE = true;

const dom = {};
[
    'canvas', 'alerts', 'inputFile', 'fileType', 'btnBrowse', 'btnUpload',
    'atList', 'selectAtomFormat', 'cdmFmtSel', 'cellDim', 'cellVec',
    'cellToMol', 'cellScale', 'stepCur', 'stepMax', 'stepSlider',
    'moleculeDropdown', 'checkboxCellEnabled',
].forEach(id => {
    dom[id] = document.getElementById(id);
});

// NOTE: `Module` needs to be declared as `var` to work with emscripten!
// noinspection ES6ConvertVarToLetConst
var Module = {
    preRun: [],
    postRun: [],
    curMol: 0,
    curStep: 0,
    canvas: setupCanvas(dom.canvas),
};

const change = {
    atoms: 1,
    cell: 2,
    fmt: 4,
    kpoints: 32,
};

change.step = change.atoms | change.cell | change.fmt;
change.mol = change.kpoints;

function setupCanvas(canvas) {
    canvas.width = canvas.clientWidth;
    canvas.height = canvas.clientHeight;

    dom.canvas.addEventListener('webglcontextlost', (e) => {
        e.preventDefault();
        if (window.confirm('Lost WebGL context. Reload page?')) {
            window.location.reload();
        }
    });

    return canvas;
}

function fillAtoms() {
    const fmt = parseInt(dom.selectAtomFormat.value);
    const at = Module.getAtomIt(Module.curMol, Module.curStep, fmt);
    const nat = Module.getNAtoms(Module.curMol, Module.curStep);
    const nbsp = '&nbsp;';

    let html = `
        <thead>
            <tr>
                <th scope="col">Type</th>
                <th scope="col">${nbsp}X</th>
                <th scope="col">${nbsp}Y</th>
                <th scope="col">${nbsp}Z</th>
            </tr>
        </thead>
    `;

    html += '<tbody>';
    for (let i = 0; i < nat; ++i) {
        const x = at.coord[0];
        const y = at.coord[1];
        const z = at.coord[2];

        html += `
            <tr data-idx=${i}>
                <td contenteditable data-idx='name' scope="row">${at.name}</td>
                <td contenteditable data-idx='0'>${x < 0 ? '' : nbsp}${x.toFixed(7)}</td>
                <td contenteditable data-idx='1'>${y < 0 ? '' : nbsp}${y.toFixed(7)}</td>
                <td contenteditable data-idx='2'>${z < 0 ? '' : nbsp}${z.toFixed(7)}</td>
            </tr>
        `;
        at.increment();
    }
    html += '</tbody>';

    dom.atList.innerHTML = html;
}

function fillCell() {
    const fmt = parseInt(dom.cdmFmtSel.value);
    const mat = Module.getCellVec(Module.curMol, Module.curStep);

    dom.cellDim.value = Module.getCellDim(Module.curMol, Module.curStep, fmt);

    for (let row = 0; row < 3; ++row) {
        for (let col = 0; col < 3; ++col) {
            dom.cellVec.rows[row + 1].cells[col + 1].innerText = mat[row][col];
        }
    }
}

function atomChanged(tgt) {
    const fmt = parseInt(dom.selectAtomFormat.value);
    const at = Module.getAtom(
        Module.curMol, Module.curStep, fmt,
        parseInt(tgt.parentElement.dataset.idx),
    );

    if (tgt.dataset.idx === 'name') {
        at.name = tgt.innerText;
    } else {
        const dir = parseInt(tgt.dataset.idx);
        const newVal = parseFloat(tgt.innerText);
        if (!Number.isNaN(newVal)) {
            const {coord} = at;
            coord[dir] = newVal;
            at.coord = coord;
        }
    }

    Module.updateView();
}

function cellDimChanged(tgt) {
    const newVal = parseFloat(tgt.value);
    if (Number.isNaN(newVal)) {
        return;
    }

    const fmt = parseInt(dom.cdmFmtSel.value);
    const scale = dom.cellScale.checked;
    const trajec = dom.cellToMol.checked;

    if (trajec) {
        for (let i = 0; i < Module.getMolNStep(Module.curMol); ++i) {
            Module.setCellDim(Module.curMol, i, newVal, fmt, scale);
        }
    } else {
        Module.setCellDim(Module.curMol, Module.curStep, newVal, fmt, scale);
    }

    Module.updateView();
    update(change.cell);
}

function cellEnabled(val) {
    const trajec = dom.cellToMol.checked;
    if (trajec) {
        for (let i = 0; i < Module.getMolNStep(Module.curMol); ++i) {
            Module.enableCell(Module.curMol, i, val);
        }
    } else {
        Module.enableCell(Module.curMol, Module.curStep, val);
    }

    Module.updateView();
}

function cellVecChanged(tgt) {
    const newVal = parseFloat(tgt.innerText);
    if (isNaN(newVal)) {
        return;
    }

    const col = tgt.dataset.idx;
    const row = tgt.parentElement.dataset.idx;
    const vec = Module.getCellVec(Module.curMol, Module.curStep);
    const trajec = document.getElementById('cellToMol').checked;
    const scale = document.getElementById('cellScale').checked;

    vec[row][col] = newVal;

    if (trajec) {
        for (let i = 0; i < Module.getMolNStep(Module.curMol); ++i) {
            Module.setCellVec(Module.curMol, i, vec, scale);
        }
    } else {
        Module.setCellVec(Module.curMol, Module.curStep, vec, scale);
    }

    Module.updateView();
    update(change.cell);
}

function readFile() {
    if (dom.inputFile.files.length !== 1) {
        return;
    }

    $('.alert').alert('close');

    const file = dom.inputFile.files[0];
    const reader = new FileReader();

    reader.readAsText(file);

    reader.onload = (e) => {
        Module.FS_createDataFile('/tmp', 'vipster.file', e.target.result, true);
        const readError = Module.readFile('/tmp/vipster.file', file.name, parseInt(dom.fileType.value));
        Module.FS_unlink('/tmp/vipster.file');

        if (readError.length) {
            $(document.body).append(createAlert('<strong>Unable to load file</strong><br>Correct format?', 'danger'));
            if (VERBOSE) {
                console.warn(readError);
            }
            return false;
        }

        // noinspection JSCheckFunctionSignatures
        const idx = parseInt($(dom.moleculeDropdown).find('a:last').data('idx')) + 1;
        const link = `<a class="dropdown-item" href="#" data-idx="${idx}">${Module.getMolName(idx)}</a>`;

        $(dom.moleculeDropdown)
            .find('.dropdown-divider')
            .show()
            .parent()
            .append(link);

        setMol(idx);
    };
    reader.readAsText(file);
}

function openFileDialogue() {
    dom.inputFile.click();
}

function setMult() {
    const x = document.getElementById('xmult').value;
    const y = document.getElementById('ymult').value;
    const z = document.getElementById('zmult').value;
    Module.setMult(parseInt(x), parseInt(y), parseInt(z));
}

function update(arg) {
    if (arg & (change.atoms | change.cell)) {
        fillAtoms();
    }

    if (arg & change.cell) {
        fillCell();
    }
}

function setStep(i) {
    Module.curStep = i;

    const hasCell = Module.hasCell(Module.curMol, Module.curStep);
    $('.if-cell').toggle(hasCell);

    dom.stepCur.innerHTML = i + 1;
    dom.selectAtomFormat.value = Module.getFmt(Module.curMol, Module.curStep);
    dom.checkboxCellEnabled.checked = hasCell;

    Module.setStep(Module.curMol, i);
    update(change.step);
}

function setMol(idx) {
    Module.curMol = idx;
    const nstep = Module.getMolNStep(idx);
    const molName = Module.getMolName(Module.curMol);

    // update molecule dropdown
    $(dom.moleculeDropdown).find('.dropdown-toggle:first').text(molName);

    dom.stepMax.innerHTML = nstep;
    dom.stepSlider.max = nstep - 1;
    dom.stepSlider.value = nstep - 1;

    $('#widget-step').toggle(nstep > 1);
    setStep(nstep - 1);
}

function createAlert(msg, type, dismissable = true) {
    return `
        <div class="alert alert-${type} ${dismissable ? 'alert-dismissible fade show' : ''}" role="alert">
            ${msg}
            <button type="button" class="close" data-dismiss="alert" aria-label="Close">
                <span aria-hidden="true">&times;</span>
            </button>
        </div>
    `;
}

$(document).ready(function () {
    // File loading
    $(dom.btnBrowse).click(openFileDialogue);
    $(dom.btnUpload).click(readFile);
    $(dom.inputFile).change(function () {
        dom.btnUpload.disabled = (this.files.length !== 1);
    });

    // Molecule loading
    $(document.body).on('click', `#${dom.moleculeDropdown.id} a`, function () {
        const idx = $(this).data('idx') || 0;
        setMol(idx);
    });

    // Hide specific cell properties if no cell present
    $(dom.checkboxCellEnabled).change(function () {
        $('.if-cell').toggle($(this).get(0).checked);
    });

    window.addEventListener('resize', function () {
        Module.canvas.width = Module.canvas.clientWidth;
        Module.canvas.height = Module.canvas.clientHeight;
    });

    $('.widget__toggle-btn').click(function () {
        $(this)
            .stop()
            .parents('.widget:first').toggleClass('closed')
            .find('.widget__body').slideToggle();
    });
});

// noinspection JSUnusedGlobalSymbols
function addParser(idx, name) {
    // eslint-disable-next-line no-undef
    $(dom.fileType).append(`<option value=${idx}>${UTF8ToString(name)}</option>`);
}
