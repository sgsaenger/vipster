#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>

#include "../common/guiwrapper.h"
#include "molecule.h"
#include "iowrapper.h"

namespace em = emscripten;
using namespace Vipster;

static GuiWrapper gui;
static std::vector<Molecule> molecules;

std::string emReadFile(std::string fn, std::string name, int fmt){
    try {
        auto d = readFile(fn, (IOFmt)fmt, name);
        molecules.push_back(d.mol);
        return "";
    } catch (std::exception &e) {
        return e.what();
    }
}
// Molecules
int emGetNMol(void){ return molecules.size();}
int emGetMolNstep(int m){ return molecules[m].getNstep();}
std::string emGetMolName(int m){ return molecules[m].getName();}

// Steps
void emSetStep(int m, int s){ gui.updateBuffers(&molecules[m].getStep(s), true); }
void emSetMult(uint8_t x, uint8_t y, uint8_t z){ gui.mult = {{x,y,z}}; }
int emGetNAtoms(int m, int s){ return molecules[m].getStep(s).getNat(); }
Atom emGetAtom(int m, int s, int fmt, int at){ return molecules[m].getStep(s).asFmt((AtomFmt)fmt)[at]; }
Step::iterator emGetAtomIt(int m, int s, int fmt){ return molecules[m].getStep(s).asFmt((AtomFmt)fmt).begin(); }
int emGetFmt(int m, int s){ return (int)molecules[m].getStep(s).getFmt();}

// Atom
std::string emGetAtName(const Atom& at){return at.name;}
void emSetAtName(Atom& at, std::string name){at.name = name;}
Vec emGetAtCoord(const Atom& at){return at.coord;}
void emSetAtCoord(Atom& at, Vec v){at.coord = v;}

// Iterator
std::string emGetItName(const Step::iterator& it){return it->name;}
void emSetItName(Step::iterator& it, std::string name){it->name = name;}
Vec emGetItCoord(const Step::iterator& it){return it->coord;}
void emSetItCoord(Step::iterator& it, Vec v){it->coord = v;}

// Cell
float emGetCellDim(int m, int s, int fmt){return molecules[m].getStep(s).getCellDim((CdmFmt)fmt);}
void emSetCellDim(int m, int s, float cdm, int fmt, bool scale){molecules[m].getStep(s).setCellDim(cdm, (CdmFmt)fmt, scale);}
Mat emGetCellVec(int m, int s){return molecules[m].getStep(s).getCellVec();}
void emSetCellVec(int m, int s, Mat vec, bool scale){molecules[m].getStep(s).setCellVec(vec, scale);}
void emEnableCell(int m, int s, bool b){molecules[m].getStep(s).enableCell(b);}
bool emHasCell(int m, int s){return molecules[m].getStep(s).hasCell();}

// Expose Canvas operations
void emUpdateView(void){ gui.updateBuffers(nullptr, true); }
void emZoom(int val){gui.zoomViewMat(val);}
void emRotate(int x, int y){gui.rotateViewMat(x,y,0);}
void emTranslate(int x, int y){gui.translateViewMat(x,y,0);}

EMSCRIPTEN_BINDINGS(vipster){
    em::function("getNMol", &emGetNMol);
    em::function("getMolNStep", &emGetMolNstep);
    em::function("getMolName", &emGetMolName);
    em::function("setStep", &emSetStep);
    em::function("setMult", &emSetMult);
    em::function("readFile", &emReadFile);
    em::function("getAtom", &emGetAtom);
    em::function("getAtomIt", &emGetAtomIt);
    em::function("getFmt", &emGetFmt);
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
    em::value_array<Vec>("Vec")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::value_array<Mat>("Mat")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::class_<Atom>("Atom")
            .property("name", &emGetAtName, &emSetAtName)
            .property("coord", &emGetAtCoord, &emSetAtCoord);
    em::class_<Step::iterator>("Step_iterator")
            .function("increment", &Step::iterator::operator++)
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
            gui.rotateViewMat(mouseEvent->canvasX-localX,
                              mouseEvent->canvasY-localY, 0);
            break;
        case OpMode::Translation:
            gui.translateViewMat(mouseEvent->canvasX-localX,
                                 -(mouseEvent->canvasY-localY), 0);
            break;
        default:
            break;
        }
        break;
    }
    localX = mouseEvent->canvasX;
    localY = mouseEvent->canvasY;
    return 1;
}

EM_BOOL wheel_event(int, const EmscriptenWheelEvent* wheelEvent, void*)
{
    gui.zoomViewMat(-wheelEvent->deltaY);
    return 1;
}

void one_iter(){
    int width, height;
    static int localWidth, localHeight;
    // handle resize
    emscripten_get_canvas_element_size("#canvas", &width, &height);
    if( width != localWidth || height != localHeight){
        gui.resizeViewMat(width, height);
        localWidth = width;
        localHeight = height;
    }
    // sync data and draw
    gui.updateVBOs();
    gui.updateViewUBO();
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
    EMSCRIPTEN_WEBGL_CONTEXT_HANDLE context = emscripten_webgl_create_context( 0, &attrs );
    if (!context)
    {
        printf("WebGL 2 is not supported!\n");
        return 0;
    }

    // init GL
    emscripten_webgl_make_context_current(context);
    glClearColor(1,1,1,1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    gui.initShaders("# version 300 es\nprecision highp float;\n", "");
    gui.initAtomVAO();
    gui.initBondVAO();
    gui.initCellVAO();
    gui.initViewUBO();
    gui.initViewMat();

    // init examples (something needs to be displayed for the renderer to not fail
    StepProper* step;
    //example H2O-vibration (crude approximation)
    molecules.emplace_back("Example Molecule", 0);
    float vibdist[] = {0,0.02,0.04,0.06,0.04,0.02,0};
    for(float f:vibdist){
        step = &molecules[0].newStep();
        step->enableCell(false);
        step->setFmt(AtomFmt::Angstrom);
        step->newAtom("H",{{(float)(-0.756+f),(float)(-0.591+f),0}});
        step->newAtom("O",{{0,0,0}});
        step->newAtom("H",{{(float)(0.756-f),(float)(-0.591+f),0}});
    }
    molecules.emplace_back("Example Crystal");
    step = &molecules[1].getStep(0);
    step->setCellDim(5.64, CdmFmt::Angstrom);
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
    int i{0};
    for(auto& fmt: IOPlugins){
        EM_ASM_({addParser($0, $1)}, i++, fmt.second->name.c_str());
    }
    EM_ASM(setMol(0));
    emscripten_set_main_loop(one_iter, 0, 1);
    gui.deleteGLObjects();
    emscripten_webgl_destroy_context(context);
    return 1;
}
