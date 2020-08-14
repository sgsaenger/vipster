Module().then(function (Module) {
    const VERBOSE = true;
    const DESKTOP_BREAKPOINT = 992;

    const dom = {};
    [
        'alerts', 'atList', 'btnDownload', 'btnUpload', 'canvas',
        'cellDim', 'cellToMol', 'cellVec', 'checkboxCellEnabled',
        'fileType', 'fileDrop', 'fileName', 'inputFile', 'loadFileModal', 'moleculeDropdown',
        'noWebGL2Modal', 'saveFileType', 'selectAtomFormat', 'stepCur', 'stepMax', 'stepSlider', 'uploadGroup'
    ].forEach(id => {
        dom[id] = document.getElementById(id);
    });

    Module.curMol = null;
    Module.curStep = null;
    Module.canvas = setupCanvas(dom.canvas);
    Module.hammer = null;
    Module.file = null;
    Module.gui = null;
    Module.molecules = [];
    
    const change = {
        atoms: 1,
        cell: 2,
        fmt: 4,
        kpoints: 8,
        param: 16,
        config: 32
    };

    change.step = change.atoms | change.cell;
    change.mol = change.kpoints;

    // initialize GUI and example molecules
    try{
        Module.gui = new Module.VipsterView('canvas');
    }catch(err){
        dom.noWebGL2Modal.style.display = "block";
    }
    Module.molecules.push(new Module.Molecule('[{"atoms":[{"coord":[-0.756,-0.591,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.756,-0.591,0.0],"name":"H"}],"fmt":"angstrom"},{"atoms":[{"coord":[-0.736,-0.571,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.736,-0.571,0.0],"name":"H"}],"fmt":"angstrom"},{"atoms":[{"coord":[-0.716,-0.5509999999999999,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.716,-0.5509999999999999,0.0],"name":"H"}],"fmt":"angstrom"},{"atoms":[{"coord":[-0.696,-0.5309999999999999,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.696,-0.5309999999999999,0.0],"name":"H"}],"fmt":"angstrom"},{"atoms":[{"coord":[-0.716,-0.5509999999999999,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.716,-0.5509999999999999,0.0],"name":"H"}],"fmt":"angstrom"},{"atoms":[{"coord":[-0.736,-0.571,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.736,-0.571,0.0],"name":"H"}],"fmt":"angstrom"},{"atoms":[{"coord":[-0.756,-0.591,0.0],"name":"H"},{"coord":[0.0,0.0,0.0],"name":"O"},{"coord":[0.756,-0.591,0.0],"name":"H"}],"fmt":"angstrom"}]'));
    Module.molecules[0].name = "Example Molecule"
    Module.molecules.push(new Module.Molecule('[{"atoms":[{"coord":[0.0,0.0,0.0],"name":"Na"},{"coord":[0.5,0.0,0.0],"name":"Cl"},{"coord":[0.5,0.5,0.0],"name":"Na"},{"coord":[0.0,0.5,0.0],"name":"Cl"},{"coord":[0.5,0.0,0.5],"name":"Na"},{"coord":[0.0,0.0,0.5],"name":"Cl"},{"coord":[0.0,0.5,0.5],"name":"Na"},{"coord":[0.5,0.5,0.5],"name":"Cl"}],"cell":{"dimension":5.64,"vectors":[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]},"fmt":"crystal"}]'));
    Module.molecules[1].name = "Example Crystal"
    setMol(0);
    for(i = 0; i < Module.nplug(); i++){
        var name = Module.plugName(i);
        if(Module.plugRead(i)){
            dom.fileType.innerHTML += `<option value=${i}>${name}</option>`
        }
        if(Module.plugWrite(i)){
            dom.saveFileType.innerHTML += `<option value=${i}>${name}</option>`
        }
    }

    function setupCanvas(canvas) {
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

    function setupHammer(canvas) {
        let hammer = new Hammer(canvas);

        hammer.get('pan').set({direction: Hammer.DIRECTION_ALL, threshold: 20});
        hammer.get('pinch').set({enable: true, threshold: 0.1});
        // hammer.get('rotate').set({enable: true});
        hammer.remove(['tap', 'doubletap', 'swipe']);

        hammer.on('press', (e) => {
            if (e.pointerType !== 'touch'){
                return;
            }

            console.log('press');
        });
        hammer.on('pressup', (e) => {
            if (e.pointerType !== 'touch'){
                return;
            }

            console.log('pressup');
        });

        hammer.on('panstart', (e) => {
            if (e.pointerType !== 'touch') {
                return;
            }

            console.log('panstart');
            this.oldX = e.deltaX;
            this.oldY = e.deltaY;
        });

        hammer.on('pan', (e) => {
            if (e.pointerType !== 'touch') {
                return;
            }

            console.log('pan');
            Module.gui.rotate(e.deltaX - this.oldX, e.deltaY - this.oldY);
            this.oldX = e.deltaX;
            this.oldY = e.deltaY;
        });

        hammer.on('pinchstart', (e) => {
            if (e.pointerType !== 'touch') {
                return;
            }

            console.log('pinchstart');
            this.oldScale = 0.0
        });
        hammer.on('pinch', (e) => {
            if (e.pointerType !== 'touch') {
                return;
            }

            console.log('pinch');
                      console.log(e.scale);
            const delta = (e.scale - this.oldScale) < 0 ? 0.98 : 1.02;
            Module.gui.zoom(delta);
            this.oldScale = e.scale;
        });

        return hammer;
    }

    function fillAtoms() {
        const at = Module.curStep.getAtomIt();
        const nat = Module.curStep.nat;
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
                    <td contenteditable data-idx='0'>${x < 0 ? '' : nbsp}${x.toFixed(6)}</td>
                    <td contenteditable data-idx='1'>${y < 0 ? '' : nbsp}${y.toFixed(6)}</td>
                    <td contenteditable data-idx='2'>${z < 0 ? '' : nbsp}${z.toFixed(6)}</td>
                </tr>
            `;
            at.increment();
        }
        html += '</tbody>';

        dom.atList.innerHTML = html;
    }

    function fillCell() {
        //const fmt = parseInt(dom.cdmFmtSel.value);
        const mat = Module.curStep.cellVec;

        dom.cellDim.value = Module.curStep.cellDim;

        for (let row = 0; row < 3; ++row) {
            for (let col = 0; col < 3; ++col) {
                dom.cellVec.rows[row + 1].cells[col + 1].innerText = mat[row][col];
            }
        }
    }

    function atomChanged(tgt) {
        const fmt = parseInt(dom.selectAtomFormat.value);
        const at = Module.curStep.getAtom( parseInt(tgt.parentElement.dataset.idx));

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

        update(change.atoms);
    }

    function cellDimChanged(tgt) {
        const newVal = parseFloat(tgt.value);
        if (Number.isNaN(newVal)) {
            return;
        }

        //const fmt = parseInt(dom.cdmFmtSel.value);
        //const scale = dom.cellScale.checked;
        const trajec = dom.cellToMol.checked;

        if (trajec) {
            for (let i = 0; i < Module.curMol.nstep; ++i) {
                Module.curMol.getStep(i).cellDim = val;
            }
        } else {
            Module.curStep.cellDim = newVal;
        }

        update(change.cell);
    }

    function cellEnabled(val) {
        const trajec = dom.cellToMol.checked;
        if (trajec) {
            for (let i = 0; i < Module.curMol.nstep; ++i) {
                Module.curMol.getStep(i).hasCell = val;
            }
        } else {
            Module.curStep.hasCell = val;
        }

        update(change.cell);
    }

    function cellVecChanged(tgt) {
        const newVal = parseFloat(tgt.innerText);
        if (isNaN(newVal)) {
            return;
        }

        const col = tgt.dataset.idx;
        const row = tgt.parentElement.dataset.idx;
        const vec = Module.curStep.cellVec;
        const trajec = dom.cellToMol.checked;
        //const scale = document.getElementById('cellScale').checked;

        vec[row][col] = newVal;

        if (trajec) {
            for (let i = 0; i < Module.curMol.nstep; ++i) {
                Module.curMol.getStep(i).cellVec = vec;
            }
        } else {
            Module.curStep.cellVec = vec;
        }

        update(change.cell);
    }

    function loadDialog() {
        dom.inputFile.value = "";
        dom.fileDrop.style.display = "";
        dom.uploadGroup.style.display = "none";
        $('#loadFileModal').modal();
    }

    function selectFile(file) {
        dom.fileDrop.style.display = "none";
        dom.uploadGroup.style.display = "";
        dom.fileName.textContent = file.name;
        dom.fileType.value = Module.guessFmt(file.name);
        Module.file = file;
    }

    function dropHandler(ev){
        ev.preventDefault();
        if(ev.dataTransfer.files.length === 1){
            selectFile(ev.dataTransfer.files[0]);
        }
    }

    function readFile() {
        $('.alert').alert('close');

        //const file = dom.inputFile.files[0];
        if(!Module.file instanceof File){
            $(document.body).append(createAlert('<strong>Trying to load something that is not a file</strong><br>Cancelling...','danger'));
            return;
        }
        const reader = new FileReader();

        reader.onload = (e) => {
            FS.createDataFile('/tmp', Module.file.name, e.target.result, true);
            try{
                Module.molecules.push(new Module.Molecule('/tmp/'+Module.file.name, parseInt(dom.fileType.value)));
            }catch(err){
                $(document.body).append(createAlert('<strong>Unable to load file</strong><br>Correct format?', 'danger'));
                if (VERBOSE) {
                    console.warn(err);
                }
                return false;
            }
            FS.unlink('/tmp/'+Module.file.name);

            // noinspection JSCheckFunctionSignatures
            const idx = Module.molecules.length - 1;
            const link = `<a class="dropdown-item" href="#" data-idx="${idx}">${Module.molecules[idx].name}</a>`;

            $(dom.moleculeDropdown)
                .find('.dropdown-divider')
                .show()
                .parent()
                .append(link);

            Module.file = null;
            setMol(idx);
        };
        reader.readAsText(Module.file);
        $('#loadFileModal').modal('hide');
    }

    function saveDialog() {
        $('#saveFileModal').modal();
    }

    function saveFile() {
        const step = parseInt(dom.stepSlider.value);
        const writeError = Module.writeFile(Module.curMol, step, parseInt(dom.saveFileType.value));
        if (writeError.length) {
            $(document.body).append(createAlert('<strong>Unable to download file</strong>', 'danger'));
            if (VERBOSE) {
                console.warn(writeError);
            }
            return false;
        }
        var data = FS.readFile('/tmp/output.file');
        FS.unlink('/tmp/output.file');
        var blob = new Blob([data.buffer], {type: "text/plain"});
        var url = window.URL.createObjectURL(blob);
        this.href = url;
        this.target = '_blank';
        this.download = Module.addExtension(Module.curMol.name, parseInt(dom.saveFileType.value));
        console.log(this);
    }

    function setMult() {
        const x = document.getElementById('xmult').value;
        const y = document.getElementById('ymult').value;
        const z = document.getElementById('zmult').value;
        Module.gui.mult = [parseInt(x), parseInt(y), parseInt(z)];
    }

    function update(arg) {
        if (arg & change.fmt) {
            Module.curStep.fmt = parseInt(dom.selectAtomFormat.value) -2;
            console.log(Module.curStep.fmt);
            if (Module.curStep.fmt < 0) {
                arg = change.cell;
            }
        }

        if (arg & (change.atoms | change.cell | change.fmt)) {
            // ensure that step-data is in a valid state
            Module.curStep.setBonds();
            // ensure that GL is in a valid state
            Module.gui.setStep(Module.curStep);
            // update atom-table
            fillAtoms();
        }

        if (arg & change.cell) {
            // update cell-widget
            fillCell();
        }
    }

    function setStep(i) {
        Module.curStep = Module.curMol.getStep(i);

        const hasCell = Module.curStep.hasCell;
        $('.if-cell').toggle(hasCell);

        dom.stepCur.innerHTML = i + 1;
        dom.selectAtomFormat.value = Module.curStep.fmt + 2;
        dom.checkboxCellEnabled.checked = hasCell;

        update(change.step);
    }

    function setMol(idx) {
        Module.curMol = Module.molecules[idx];
        const nstep = Module.curMol.nstep;
        const molName = Module.curMol.name;

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

    function resizeCanvas() {
        Module.canvas.width = Module.canvas.clientWidth;
        Module.canvas.height = Module.canvas.clientHeight;

        // Might be hidden right now after resizing back into desktop viewport
        // because of mobile menu logic
        if ($(window).width() >= DESKTOP_BREAKPOINT) {
            $('main').show();
        }
    }

    $(document).ready(function () {
        Module.hammer = setupHammer(canvas);

        // Set correct canvas size on resize
        window.addEventListener('resize', resizeCanvas);

        dom.btnDownload.onclick = saveFile;

        // Molecule selecting
        $(document.body).on('click', `#${dom.moleculeDropdown.id} a`, function () {
            const idx = $(this).data('idx') || 0;
            setMol(idx);
        });

        // Hide specific cell properties if no cell present
        $(dom.checkboxCellEnabled).change(function () {
            $('.if-cell').toggle($(this).get(0).checked);
        });

        $('.widget_toggle').click(function () {
            $(this)
                .stop()
                .parents('.widget:first').toggleClass('closed')
                .find('.widget__body').slideToggle();
        });

        const main = $('main');
        $('#controls__collapse')
            .on('show.bs.collapse', () => main.hide())
            .on('hide.bs.collapse', () => {
                main.show();
                resizeCanvas();
            });
    });
});
