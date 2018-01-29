var Module = {
    preRun: [],
    postRun: [],
    curMol: 0,
    curStep: 0,
    canvas: (function() {
      var canvas = document.getElementById('canvas');
      canvas.addEventListener("webglcontextlost", function(e) { alert('WebGL context lost. You will need to reload the page.'); e.preventDefault(); }, false);
      canvas.width = canvas.clientWidth;
      canvas.height = canvas.clientHeight;
      return canvas;
    })(),
}

function fillAtoms(){
    var at = Module.getAtoms(Module.curMol, Module.curStep);
    var nat = Module.getNAtoms(Module.curMol, Module.curStep);
    var atList = document.getElementById('atList')
    atList.innerHTML = "<tr><th>Type</th><th>X</th><th>Y</th><th>Z</th></tr>";
    for(i=0; i<nat; ++i){
        atList.innerHTML += "<tr> \
<td contenteditable onInput=nameChanged("+i+",this.innerText)>" + at.name + "</td> \
<td contenteditable onInput=coordChanged("+i+",this.parentElement)>" + at.coord[0] + "</td> \
<td contenteditable onInput=coordChanged("+i+",this.parentElement)>" + at.coord[1] + "</td> \
<td contenteditable onInput=coordChanged("+i+",this.parentElement)>" + at.coord[2] + "</td></tr>";
        at.increment();
    }
}

function nameChanged(idx, name){
    Module.getAtom(Module.curMol, Module.curStep, idx).name = name;
}

function coordChanged(idx, line){
    var x = parseFloat(line.children[1].innerText);
    var y = parseFloat(line.children[2].innerText);
    var z = parseFloat(line.children[3].innerText);
    if(!isNaN(x)&&!isNaN(y)&&!isNaN(z)){
        Module.getAtom(Module.curMol, Module.curStep, idx).coord = [x,y,z];
    }
}

function readFile(e){
    var file = document.getElementById('upfile').files[0];
    var reader = new FileReader();
    var molList = document.getElementById('molList');
    reader.onload = function(e){
        molList.innerHTML += ('<li onclick="setMol(getMolListIdx(event.target))">'+file.name+"</li>");
        Module.FS_createDataFile("/tmp","test.file", e.target.result, true);
        Module.readFile("/tmp/test.file", file.name, parseInt(document.getElementById('uptype').value));
        Module.FS_unlink("/tmp/test.file")
        setMol(molList.childElementCount-1);
    };
    reader.readAsText(file);
    document.getElementById('output').innerHTML = "Success!";
}

function setMult(e){
    var x = document.getElementById('xmult').value;
    var y = document.getElementById('ymult').value;
    var z = document.getElementById('zmult').value;
    Module.setMult(parseInt(x),parseInt(y),parseInt(z));
}

function getMolListIdx(li){
    var molList = li.parentElement.children;
    for(i=0;i<li.parentElement.childElementCount;i++){
        if (li == molList[i]) break;
    }
    return i;
}

function setMol(i){
    Module.curMol = i;
    var nstep = Module.getMolNStep(i);
    var slider = document.getElementById('stepSlider');
    document.getElementById('stepMax').innerHTML = nstep;
    document.getElementById('stepCur').innerHTML = 1;
    slider.value = 0;
    slider.max = nstep-1;
    Module.curStep = 0;
    Module.setStep(i, 0);
    fillAtoms();
}

function setStep(e){
    var i = parseInt(e.target.value);
    document.getElementById('stepCur').innerHTML = i+1;
    Module.curStep = i;
    Module.setStep(Module.curMol, i);
    fillAtoms();
}
