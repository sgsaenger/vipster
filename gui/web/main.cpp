#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>

#include "guiwrapper.h"
#include "vipster/molecule.h"
#include "vipster/fileio.h"

namespace em = emscripten;
using namespace Vipster;

static GUI::GlobalData glGlobals{};
static GuiWrapper gui{glGlobals, settings};
static std::vector<Molecule> molecules;

static PluginList plugins = defaultPlugins();

std::string emReadFile(std::string fn, int fmt){
    try {
        auto [m, p, d] = readFile(fn, plugins[fmt]);
        molecules.push_back(std::move(m));
        return "";
    } catch (std::exception &e) {
        return e.what();
    }
}

std::string emWriteFile(int m, int s, int f){
    try{
        const auto& plug = plugins[f];
        Parameter param{nullptr};
        if(plug->makeParam){
            param = plug->makeParam();
        }
        Preset preset{nullptr};
        if(plug->makePreset){
            preset = plug->makePreset();
        }
        writeFile("/tmp/output.file", plug, molecules[m],
                  (size_t)s, param, preset);
        return "";
    } catch(std::exception &e) {
        return e.what();
    }
}

// Molecules
int emGetNMol(void){ return molecules.size();}
int emGetMolNstep(int m){ return molecules[m].getNstep();}
std::string emGetMolName(int m){ return molecules[m].name;}

// Steps
void emSetStep(int m, int s){ gui.setMainStep(&molecules[m].getStep(s)); }
void emSetMult(uint8_t x, uint8_t y, uint8_t z){ gui.mult = {{x,y,z}}; }
int emGetNAtoms(){ return gui.curStep->getNat(); }
Step::atom emGetAtom(int at){ return (*gui.curStep)[at]; }
Step::iterator emGetAtomIt(){ return gui.curStep->begin(); }
int emGetFmt(){ return (int)gui.curStep->getFmt(); }
void emSetFmt(int fmt){ gui.curStep->setFmt((AtomFmt) fmt, true); }

// Atom
std::string emGetAtName(const Step::atom& at){return at.name;}
void emSetAtName(Step::atom& at, const std::string &name){at.name = name;}
Vec emGetAtCoord(const Step::atom& at){return at.coord;}
void emSetAtCoord(Step::atom& at, Vec v){at.coord = v;}

// Iterator
std::string emGetItName(const Step::iterator& it){return it->name;}
void emSetItName(Step::iterator& it, const std::string &name){it->name = name;}
Vec emGetItCoord(const Step::iterator& it){return it->coord;}
void emSetItCoord(Step::iterator& it, Vec v){it->coord = v;}

// Cell
double emGetCellDim(int fmt){return gui.curStep->getCellDim((AtomFmt)fmt);}
void emSetCellDim(double cdm, int fmt, bool scale){gui.curStep->setCellDim(cdm, (AtomFmt)fmt, scale);}
Mat emGetCellVec(){return gui.curStep->getCellVec();}
void emSetCellVec(Mat vec, bool scale){ gui.curStep->setCellVec(vec, scale);}
void emEnableCell(bool b){ gui.curStep->enableCell(b); }
bool emHasCell(){return gui.curStep->hasCell(); }

// Expose Canvas operations
void emUpdateView(void){ gui.updateMainStep(); }
void emZoom(float val){gui.zoomViewMat(val);}
void emRotate(int x, int y){gui.rotateViewMat(x,y,0);}
void emTranslate(int x, int y){gui.translateViewMat(x,y,0);}

// validate Cache
void emEvalBonds(void){ gui.curStep->setBonds(); }

int emGuessFmt(std::string file){
    auto pos = file.find_last_of('.');
    if(pos != file.npos){
        std::string ext = file.substr(pos+1);
        auto pos = std::find_if(plugins.begin(), plugins.end(),
                                [&](const Plugin* plug){
                                    return plug->extension == ext;
                                });
        if(pos == plugins.end()){
            return 0;
        }else{
            return std::distance(plugins.begin(), pos);
        }
    }else{
        return 0;
    }
};

std::string emFmtName(int m, int f){
    auto name = molecules[m].name;
    auto dot = name.find_last_of('.');
    if(dot != name.npos){
        name = name.substr(0, dot);
    }
    return name + '.' + plugins[f]->extension;
}


EMSCRIPTEN_BINDINGS(vipster){
    em::function("evalBonds", &emEvalBonds);
    em::function("getNMol", &emGetNMol);
    em::function("getMolNStep", &emGetMolNstep);
    em::function("getMolName", &emGetMolName);
    em::function("setStep", &emSetStep);
    em::function("setMult", &emSetMult);
    em::function("readFile", &emReadFile);
    em::function("writeFile", &emWriteFile);
    em::function("getAtom", &emGetAtom);
    em::function("getAtomIt", &emGetAtomIt);
    em::function("getFmt", &emGetFmt);
    em::function("setFmt", &emSetFmt);
    em::function("getNAtoms", &emGetNAtoms);
    em::function("getCellDim", &emGetCellDim);
    em::function("setCellDim", &emSetCellDim);
    em::function("getCellVec", &emGetCellVec);
    em::function("setCellVec", &emSetCellVec);
    em::function("enableCell", &emEnableCell);
    em::function("hasCell", &emHasCell);
    em::function("updateView", &emUpdateView);
    em::function("zoom", &emZoom);
    em::function("rotate", &emRotate);
    em::function("translate", &emTranslate);
    em::function("guessFmt", &emGuessFmt);
    em::function("getFormattedName", &emFmtName);
    em::value_array<Vec>("Vec")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::value_array<Mat>("Mat")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::class_<Step::atom>("Atom")
            .property("name", &emGetAtName, &emSetAtName)
            .property("coord", &emGetAtCoord, &emSetAtCoord);
    em::class_<Step::iterator>("Step_iterator")
            .function("increment", em::select_overload<Step::iterator&()>(&Step::iterator::operator++))
            .property("name", &emGetItName, &emSetItName)
            .property("coord", &emGetItCoord, &emSetItCoord);
}

EM_BOOL mouse_event(int eventType, const EmscriptenMouseEvent* mouseEvent, void*)
{
    enum class MouseMode { Camera, Select, Modify};
    enum class OpMode { None, Rotation, Translation };
    static OpMode currentOp = OpMode::None;
    static long localX, localY;
    switch (eventType) {
    case EMSCRIPTEN_EVENT_MOUSEDOWN:
        if(currentOp == OpMode::None){
            int button = mouseEvent->button | mouseEvent->altKey | mouseEvent->ctrlKey << 1;
            switch(button){
            case 0:
                currentOp = OpMode::Rotation;
                break;
            case 1:
                currentOp = OpMode::Translation;
                break;
            case 2:
                gui.alignViewMat(GuiWrapper::alignDir::z);
            }
        }
        break;
    case EMSCRIPTEN_EVENT_MOUSEUP:
        if(currentOp!=OpMode::None && !mouseEvent->buttons){
            currentOp = OpMode::None;
        }
        break;
    case EMSCRIPTEN_EVENT_MOUSELEAVE:
        currentOp = OpMode::None;
        break;
    case EMSCRIPTEN_EVENT_MOUSEMOVE:
        switch(currentOp){
        case OpMode::Rotation:
            gui.rotateViewMat(mouseEvent->clientX-localX,
                              mouseEvent->clientY-localY, 0);
            break;
        case OpMode::Translation:
            gui.translateViewMat(mouseEvent->clientX-localX,
                                 -(mouseEvent->clientY-localY), 0);
            break;
        default:
            break;
        }
        break;
    }
    localX = mouseEvent->clientX;
    localY = mouseEvent->clientY;
    return 1;
}

EM_BOOL wheel_event(int, const EmscriptenWheelEvent* wheelEvent, void*)
{
    gui.zoomViewMat(wheelEvent->deltaY<0?1.1:0.9);
    return 1;
}

void main_loop(){
    int width, height;
    static int localWidth, localHeight;
    // handle resize
    emscripten_get_canvas_element_size("#canvas", &width, &height);
    if( width != localWidth || height != localHeight){
        gui.resizeViewMat(width, height);
        localWidth = width;
        localHeight = height;
    }
    // draw
    gui.draw();
}

int main()
{
    // create WebGL2 context
    EmscriptenWebGLContextAttributes attrs;
    emscripten_webgl_init_context_attributes(&attrs);
    attrs.enableExtensionsByDefault = 1;
    attrs.majorVersion = 2;
    attrs.minorVersion = 0;
    EMSCRIPTEN_WEBGL_CONTEXT_HANDLE context = emscripten_webgl_create_context( "#canvas", &attrs );
    if (!context)
    {
        printf("WebGL 2 is not supported!\n");
        EM_ASM(alertWebGL());
        return 0;
    }

    // init GL
    emscripten_webgl_make_context_current(context);
    gui.initGL();

    // init examples (something needs to be displayed for the renderer to not fail
    Step* step;
    //example H2O-vibration (crude approximation)
    molecules.emplace_back("Example Molecule", 0);
    double vibdist[] = {0,0.02,0.04,0.06,0.04,0.02,0};
    for(double f:vibdist){
        step = &molecules[0].newStep();
        step->enableCell(false);
        step->setFmt(AtomFmt::Angstrom);
        step->newAtom("H",{{-0.756+f,-0.591+f,0}});
        step->newAtom("O",{{0,0,0}});
        step->newAtom("H",{{0.756-f,-0.591+f,0}});
    }
    molecules.emplace_back("Example Crystal");
    step = &molecules[1].getStep(0);
    step->setCellDim(5.64, AtomFmt::Angstrom);
    step->setFmt(AtomFmt::Crystal);
    step->newAtom("Na",{{0.0,0.0,0.0}});
    step->newAtom("Cl",{{0.5,0.0,0.0}});
    step->newAtom("Na",{{0.5,0.5,0.0}});
    step->newAtom("Cl",{{0.0,0.5,0.0}});
    step->newAtom("Na",{{0.5,0.0,0.5}});
    step->newAtom("Cl",{{0.0,0.0,0.5}});
    step->newAtom("Na",{{0.0,0.5,0.5}});
    step->newAtom("Cl",{{0.5,0.5,0.5}});

    // handle input
    emscripten_set_mousedown_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_mouseup_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_mousemove_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_mouseleave_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_wheel_callback("#canvas", nullptr, 0, wheel_event);

    //start
    for(int i=0; i<plugins.size(); ++i){
        EM_ASM_({addParser($0, $1)}, i, plugins[i]->name.c_str());
        if(plugins[i]->writer){
            EM_ASM_({addWriter($0, $1)}, i, plugins[i]->name.c_str());
        }
    }
    EM_ASM(setMol(0));
    emscripten_set_main_loop(main_loop, 0, 1);
    emscripten_webgl_destroy_context(context);
    return 1;
}
