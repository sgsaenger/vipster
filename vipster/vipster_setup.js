const dom = {};
[ 'canvas', 'atList', 'atFmtSel', 'cdmFmtSel', 'cellDim', 'cellVec',
  'cellToMol', 'cellScale', 'stepCur', 'stepMax', 'stepSlider',
].forEach(id => {
  dom[id] = document.getElementById(id);
});

// NOTE: this needs to be declared as `var` to work with emscripten!
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
  const fmt = parseInt(dom.atFmtSel.value);
  const at = Module.getAtomIt(Module.curMol, Module.curStep, fmt);
  const nat = Module.getNAtoms(Module.curMol, Module.curStep);

  let html = `
    <thead>
      <tr>
        <th scope="col">Type</th>
        <th scope="col">X</th>
        <th scope="col">Y</th>
        <th scope="col">Z</th>
      </tr>
    </thead>
  `;

  html += '<tbody>';
  for (let i = 0; i < nat; ++i) {
    html += `
      <tr data-idx=${i}>
        <td contenteditable data-idx='name' scope="row">${at.name}</td>
        <td contenteditable data-idx='0'>${at.coord[0]}</td>
        <td contenteditable data-idx='1'>${at.coord[1]}</td>
        <td contenteditable data-idx='2'>${at.coord[2]}</td>
      </tr>
    `;
    at.increment();
  }
  html += '</tbody>';

  dom.atList.innerHTML = html;
}

function fillCell() {
  console.log(dom.cdmFmtSel);
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
  const fmt = parseInt(dom.atFmtSel.value);
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
  vec[row][col] = newVal;
  const trajec = document.getElementById('cellToMol').checked;
  const scale = document.getElementById('cellScale').checked;
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
  const file = document.getElementById('upfile').files[0];
  const molList = document.getElementById('molList');
  const reader = new FileReader();
  reader.onload = (e) => {
    Module.FS_createDataFile('/tmp', 'test.file', e.target.result, true);
    Module.readFile(
      '/tmp/test.file', file.name,
      Number.parseInt(document.getElementById('uptype').value),
    );
    Module.FS_unlink('/tmp/test.file');
    molList.innerHTML += `<li onclick="setMol(getMolListIdx(event.target))">${file.name}</li>`;
    setMol(molList.childElementCount - 1);
  };
  reader.readAsText(file);
}

function setMult() {
  const x = document.getElementById('xmult').value;
  const y = document.getElementById('ymult').value;
  const z = document.getElementById('zmult').value;
  Module.setMult(parseInt(x), parseInt(y), parseInt(z));
}

function update(change) {
  if (change & (change.atoms | change.cell)) {
    fillAtoms();
  }

  if (change & change.cell) {
    fillCell();
  }
}

function setStep(i) {
  dom.stepCur.innerHTML = i + 1;
  Module.curStep = i;
  Module.setStep(Module.curMol, i);
  dom.atFmtSel.value = Module.getFmt(Module.curMol, Module.curStep);
  update(change.step);
}

function setMol(idx) {
  Module.curMol = idx;
  const nstep = Module.getMolNStep(idx);
  dom.stepMax.innerHTML = nstep;
  dom.stepSlider.max = nstep - 1;
  dom.stepSlider.value = nstep - 1;
  setStep(nstep - 1);
}

$(document).ready(function () {
  const widgets = {
    moleculeList: $('.widget-molecule-list:first'),
  };

  widgets.moleculeList.find('a').click(function () {
    const idx = $(this).data('idx') || 0;
    setMol(idx);
  });
});
