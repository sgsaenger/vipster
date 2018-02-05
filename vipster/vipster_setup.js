const change = {
  atoms: 1,
  cell: 2,
  fmt: 4,
  kpoints: 32,
};
change.step = change.atoms | change.cell | change.fmt;
change.mol = change.kpoints;

function setupCanvas() {
  const canvas = document.getElementById('canvas');

  canvas.width = canvas.clientWidth;
  canvas.height = canvas.clientHeight;

  canvas.addEventListener('webglcontextlost', (e) => {
    e.preventDefault();
    if (window.confirm('Lost WebGL context. Reload page?')) {
      window.location.reload();
    }
  });

  return canvas;
}

var Module = {
  preRun: [],
  postRun: [],
  curMol: 0,
  curStep: 0,
  canvas: setupCanvas(),
};

function fillAtoms() {
  const fmt = parseInt(document.getElementById('atFmtSel').value, 10);
  const at = Module.getAtomIt(Module.curMol, Module.curStep, fmt);
  const nat = Module.getNAtoms(Module.curMol, Module.curStep);
  const atList = document.getElementById('atList');

  let html = '<tr><th>Type</th><th>X</th><th>Y</th><th>Z</th></tr>';

  for (let i = 0; i < nat; ++i) {
    html += `
      <tr data-idx=${i}>
        <td contenteditable data-idx='name'>${at.name}</td>
        <td contenteditable data-idx='0'>${at.coord[0]}</td>
        <td contenteditable data-idx='1'>${at.coord[1]}</td>
        <td contenteditable data-idx='2'>${at.coord[2]}</td>
      </tr>
    `;

    atList.innerHTML = html;
    at.increment();
  }
}

function fillCell() {
  const cdm = document.getElementById('cellDim');
  const fmt = parseInt(document.getElementById('cdmFmtSel').value, 10);
  cdm.value = Module.getCellDim(Module.curMol, Module.curStep, fmt);
  const mat = Module.getCellVec(Module.curMol, Module.curStep);
  const vec = document.getElementById('cellVec');
  for (let row = 0; row < 3; ++row) {
    for (let col = 0; col < 3; ++col) {
      vec.rows[row + 1].cells[col + 1].innerText = mat[row][col];
    }
  }
}

function atomChanged(tgt) {
  const fmt = parseInt(document.getElementById('atFmtSel').value, 10);
  const at = Module.getAtom(
    Module.curMol, Module.curStep, fmt,
    parseInt(tgt.parentElement.dataset.idx, 10),
  );
  if (tgt.dataset.idx === 'name') {
    at.name = tgt.innerText;
  } else {
    const dir = parseInt(tgt.dataset.idx, 10);
    const newVal = parseFloat(tgt.innerText);
    if (!Number.isNaN(newVal)) {
      const { coord } = at;
      coord[dir] = newVal;
      at.coord = coord;
    }
  }
  Module.updateView();
}

function cellDimChanged(tgt) {
  const newVal = parseFloat(tgt.value);
  if (Number.isNaN(newVal)) { return; }
  const fmt = parseInt(document.getElementById('cdmFmtSel').value, 10);
  const scale = document.getElementById('cellScale').checked;
  const trajec = document.getElementById('cellToMol').checked;
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
  const trajec = document.getElementById('cellToMol').checked;
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
  if (isNaN(newVal)) { return; }
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
      Number.parseInt(document.getElementById('uptype').value, 10),
    );
    Module.FS_unlink('/tmp/test.file');
    molList.innerHTML += `<li onclick="setMol(getMolListIdx(event.target))">${file.name}</li>`;
    setMol(molList.childElementCount - 1);
  };
  reader.readAsText(file);
}

function setMult(e) {
  const x = document.getElementById('xmult').value;
  const y = document.getElementById('ymult').value;
  const z = document.getElementById('zmult').value;
  Module.setMult(parseInt(x), parseInt(y), parseInt(z));
}

function getMolListIdx(li) {
  const molList = li.parentElement.children;
  let i;
  for (i = 0; i < li.parentElement.childElementCount; i++) {
    if (li === molList[i]) { break; }
  }
  return i;
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
  document.getElementById('stepCur').innerHTML = i + 1;
  Module.curStep = i;
  Module.setStep(Module.curMol, i);
  document.getElementById('atFmtSel').value = Module.getFmt(Module.curMol, Module.curStep);
  update(change.step);
}

function setMol(i) {
  Module.curMol = i;
  const nstep = Module.getMolNStep(i);
  const slider = document.getElementById('stepSlider');
  document.getElementById('stepMax').innerHTML = nstep;
  slider.max = nstep - 1;
  slider.value = nstep - 1;
  setStep(nstep - 1);
}
