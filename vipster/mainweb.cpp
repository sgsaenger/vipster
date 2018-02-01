#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>

#include "guiwrapper.h"
#include "molecule.h"
#include "atomproper.h"
#include "iowrapper.h"
#include "atom_model.h"

namespace em = emscripten;
using namespace Vipster;

static GuiWrapper gui;
static std::vector<Molecule> molecules;

void emReadFile(std::string fn, std::string name, int fmt){
    auto d = readFile(fn, (IOFmt)fmt, name);
    molecules.push_back(d->mol);
}
// Molecules
int emGetNMol(void){ return molecules.size();}
int emGetMolNstep(int m){ return molecules[m].getNstep();}
std::string emGetMolName(int m){ return molecules[m].getName();}
// Steps
void emSetStep(int m, int s){ gui.updateBuffers(&molecules[m].getStep(s), true); }
void emUpdateView(void){ gui.updateBuffers(nullptr, true); }
void emSetMult(uint8_t x, uint8_t y, uint8_t z){ gui.mult = {{x,y,z}}; }
int emGetNAtoms(int m, int s){ return molecules[m].getStep(s).getNat(); }
AtomRef emGetAtom(int m, int s, int fmt, int at){ return molecules[m].getStep(s).asFmt((AtomFmt)fmt)[at]; }
Step::iterator emGetAtomIt(int m, int s, int fmt){ return molecules[m].getStep(s).asFmt((AtomFmt)fmt).begin(); }
// Atom
std::string emGetAtName(const AtomRef& at){return at.name;}
void emSetAtName(AtomRef& at, std::string name){at.name = name;}
Vec emGetAtCoord(const AtomRef& at){return at.coord;}
void emSetAtCoord(AtomRef& at, Vec v){at.coord = v;}
// Iterator
std::string emGetItName(const Step::iterator& it){return (*it).name;}
void emSetItName(Step::iterator& it, std::string name){(*it).name = name;}
Vec emGetItCoord(const Step::iterator& it){return (*it).coord;}
void emSetItCoord(Step::iterator& it, Vec v){(*it).coord = v;}
// Cell
float emGetCellDim(int m, int s, int fmt){return molecules[m].getStep(s).getCellDim((CdmFmt)fmt);}
void emSetCellDim(int m, int s, float cdm, int fmt, bool scale){molecules[m].getStep(s).setCellDim(cdm, (CdmFmt)fmt, scale);}
Mat emGetCellVec(int m, int s){return molecules[m].getStep(s).getCellVec();}
void emSetCellVec(int m, int s, Mat vec, bool scale){molecules[m].getStep(s).setCellVec(vec, scale);}
void emEnableCell(int m, int s, bool b){molecules[m].getStep(s).enableCell(b);}

EMSCRIPTEN_BINDINGS(vipster){
    em::function("getNMol", &emGetNMol);
    em::function("getMolNStep", &emGetMolNstep);
    em::function("getMolName", &emGetMolName);
    em::function("setStep", &emSetStep);
    em::function("setMult", &emSetMult);
    em::function("readFile", &emReadFile);
    em::function("getAtom", &emGetAtom);
    em::function("getAtomIt", &emGetAtomIt);
    em::function("getNAtoms", &emGetNAtoms);
    em::function("getCellDim", &emGetCellDim);
    em::function("setCellDim", &emSetCellDim);
    em::function("getCellVec", &emGetCellVec);
    em::function("setCellVec", &emSetCellVec);
    em::function("enableCell", &emEnableCell);
    em::function("updateView", &emUpdateView);
    em::value_array<Vec>("Vec")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::value_array<Mat>("Mat")
            .element(em::index<0>())
            .element(em::index<1>())
            .element(em::index<2>());
    em::class_<AtomRef>("AtomRef")
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

EM_BOOL touch_event(int eventType, const EmscriptenTouchEvent* touchEvent, void*)
{
    constexpr long translateDelta = 10, scaleDelta = 10;
    enum class TouchMode { None, Rotate, Scale, Translate};
    static TouchMode tMode = TouchMode::None;
    static long local1X, local2X, local1Y, local2Y, distance, transX, transY;
    long tmp=0, tmp2=0;
    switch (eventType) {
    case EMSCRIPTEN_EVENT_TOUCHSTART:
        if(tMode != TouchMode::None) break;
        switch(touchEvent->numTouches){
        case 2:
            local2X = touchEvent->touches[1].canvasX;
            local2Y = touchEvent->touches[1].canvasY;
            distance = std::sqrt(std::pow(local2X-touchEvent->touches[0].canvasX,2)+
                                 std::pow(local2Y-touchEvent->touches[0].canvasY,2));
            transX = (local2X + touchEvent->touches[0].canvasX)/2;
            transY = (local2Y + touchEvent->touches[0].canvasY)/2;
        case 1:
            tMode = TouchMode::Rotate;
            local1X = touchEvent->touches[0].canvasX;
            local1Y = touchEvent->touches[0].canvasY;
            break;
        default:
            gui.alignViewMat(GuiWrapper::alignDir::z);
            break;
        }
        break;
    case EMSCRIPTEN_EVENT_TOUCHMOVE:
        switch(touchEvent->numTouches){
        case 1:
            if(tMode != TouchMode::Rotate) break;
            gui.rotateViewMat(touchEvent->touches[0].canvasX-local1X,
                    touchEvent->touches[0].canvasY-local1Y,0);
            local1X = touchEvent->touches[0].canvasX;
            local1Y = touchEvent->touches[0].canvasY;
            break;
        case 2:
            if(tMode < TouchMode::Scale){
                tmp = std::sqrt(std::pow(touchEvent->touches[1].canvasX-
                                         touchEvent->touches[0].canvasX,2)+
                                std::pow(touchEvent->touches[1].canvasY-
                                         touchEvent->touches[0].canvasY,2));
                if (std::abs(tmp-distance)>scaleDelta){
                    tMode = TouchMode::Scale;
                }else if(local1X-touchEvent->touches[0].canvasX > translateDelta||
                         local1Y-touchEvent->touches[0].canvasY > translateDelta||
                         local2X-touchEvent->touches[1].canvasX > translateDelta||
                         local2Y-touchEvent->touches[1].canvasY > translateDelta){
                    tMode = TouchMode::Translate;
                }
            }
            if(tMode < TouchMode::Scale) break;
            if(tMode == TouchMode::Scale){
                if(!tmp)
                    tmp = std::sqrt(std::pow(touchEvent->touches[1].canvasX-
                                             touchEvent->touches[0].canvasX,2)+
                                    std::pow(touchEvent->touches[1].canvasY-
                                             touchEvent->touches[0].canvasY,2));
                gui.zoomViewMat(tmp-distance);
                distance = tmp;
            }else if(tMode == TouchMode::Translate){
                tmp = (touchEvent->touches[0].canvasX + touchEvent->touches[1].canvasX)/2;
                tmp2 = (touchEvent->touches[0].canvasY + touchEvent->touches[1].canvasY)/2;
                gui.translateViewMat(tmp-transX, -(tmp2-transY), 0);
                transX = tmp;
                transY = tmp2;
            }
            local1X = touchEvent->touches[0].canvasX;
            local1Y = touchEvent->touches[0].canvasY;
            local2X = touchEvent->touches[1].canvasX;
            local2Y = touchEvent->touches[1].canvasY;
            break;
        }
        break;
    case EMSCRIPTEN_EVENT_TOUCHEND:
    case EMSCRIPTEN_EVENT_TOUCHCANCEL:
        tMode = TouchMode::None;
        break;
    }
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
        step->newAtom(AtomProper{"H",{{(float)(-0.756+f),(float)(-0.591+f),0}}});
        step->newAtom(AtomProper{"O",{{0,0,0}}});
        step->newAtom(AtomProper{"H",{{(float)(0.756-f),(float)(-0.591+f),0}}});
    }
    molecules.emplace_back("Example Crystal");
    step = &molecules[1].getStep(0);
    step->setCellDim(5.64, CdmFmt::Angstrom);
    step->setFmt(AtomFmt::Crystal);
    step->newAtom(AtomProper{"Na",{{0.0,0.0,0.0}}});
    step->newAtom(AtomProper{"Cl",{{0.5,0.0,0.0}}});
    step->newAtom(AtomProper{"Na",{{0.5,0.5,0.0}}});
    step->newAtom(AtomProper{"Cl",{{0.0,0.5,0.0}}});
    step->newAtom(AtomProper{"Na",{{0.5,0.0,0.5}}});
    step->newAtom(AtomProper{"Cl",{{0.0,0.0,0.5}}});
    step->newAtom(AtomProper{"Na",{{0.0,0.5,0.5}}});
    step->newAtom(AtomProper{"Cl",{{0.5,0.5,0.5}}});

    // handle input
    emscripten_set_mousedown_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_mouseup_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_mousemove_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_mouseleave_callback("#canvas", nullptr, 0, mouse_event);
    emscripten_set_wheel_callback("#canvas", nullptr, 0, wheel_event);
    emscripten_set_touchstart_callback("#canvas", nullptr, 0, touch_event);
    emscripten_set_touchmove_callback("#canvas", nullptr, 0, touch_event);
    emscripten_set_touchend_callback("#canvas", nullptr, 0, touch_event);

    //start
    EM_ASM(setMol(0));
    emscripten_set_main_loop(one_iter, 0, 1);
    gui.deleteGLObjects();
    emscripten_webgl_destroy_context(context);
    return 1;
}
