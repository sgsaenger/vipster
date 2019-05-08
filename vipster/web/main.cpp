#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>
#include <emscripten/vr.h>

#include "../common/guiwrapper.h"
#include "configfile.h"
#include "molecule.h"
#include "io.h"

namespace em = emscripten;
using namespace Vipster;

static GuiWrapper gui;
static std::vector<Molecule> molecules;

static VRDisplayHandle handle;
static unsigned long vrWidth, vrHeight;
static bool vrMoving{false}, vrHasPos{false};
static Vec vrPos{0,0,-10};

std::string emReadFile(std::string fn, std::string name, int fmt){
    try {
        auto d = readFile(fn, (IOFmt)fmt, name);
        molecules.push_back(d.mol);
        return "";
    } catch (std::exception &e) {
        return e.what();
    }
}

std::string emWriteFile(int m, int s, int f){
    try{
        auto fmt = static_cast<IOFmt>(f);
        const auto& plug = IOPlugins.at(fmt);
        std::unique_ptr<IO::BaseParam> param{nullptr};
        if(plug->arguments & IO::Plugin::Param){
            param = plug->makeParam("");
        }
        std::unique_ptr<IO::BaseConfig> config{nullptr};
        if(plug->arguments & IO::Plugin::Config){
            config = plug->makeConfig("");
        }
        writeFile("/tmp/output.file", fmt, molecules[m],
                  param.get(), config.get(), {(size_t)s});
        return "";
    } catch(std::exception &e) {
        return e.what();
    }
}

// Molecules
int emGetNMol(void){ return molecules.size();}
int emGetMolNstep(int m){ return molecules[m].getNstep();}
std::string emGetMolName(int m){ return molecules[m].getName();}

// Steps
void emSetStep(int m, int s){ gui.setMainStep(&molecules[m].getStep(s), true); }
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
void emUpdateView(void){ gui.updateMainStep(true); }
void emZoom(float val){gui.zoomViewMat(val);}
void emRotate(int x, int y){gui.rotateViewMat(x,y,0);}
void emTranslate(int x, int y){gui.translateViewMat(x,y,0);}

// validate Cache
void emEvalCache(void){ gui.curStep->evaluateCache(); }

int emGuessFmt(std::string file){
    auto pos = file.find_last_of('.');
    if(pos != file.npos){
        std::string ext = file.substr(pos+1);
        auto pos = std::find_if(IOPlugins.begin(), IOPlugins.end(),
                                [&](const decltype(IOPlugins)::value_type& pair){
                                    return pair.second->extension == ext;
                                });
        if(pos == IOPlugins.end()){
            return 0;
        }else{
            return static_cast<int>(pos->first);
        }
    }else{
        return 0;
    }
};

std::string emFmtName(int m, int f){
    auto name = molecules[m].getName();
    auto dot = name.find_last_of('.');
    if(dot != name.npos){
        name = name.substr(0, dot);
    }
    return name + '.' + IOPlugins.at((IOFmt)f)->extension;
}

void emVrMove(int val){vrMoving = val;}

EMSCRIPTEN_BINDINGS(vipster){
    em::function("evalCache", &emEvalCache);
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
    em::function("vrToggleMove", &emVrMove);
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

void vr_update_pos(const VRFrameData& data){
//    if(vrHasPos){
//        vrPos[0] = data.pose.position.x;
//        vrPos[1] = data.pose.position.y;
//        vrPos[2] = data.pose.position.z;
//    }else if(vrMoving){
    if(vrMoving && !vrHasPos){
        Vec tmp{
            data.pose.orientation.x,
            data.pose.orientation.y,
            data.pose.orientation.z,
        };
        tmp /= Vec_length(tmp);
        vrPos += tmp;
    }
}

void vr_loop(){
    if (!emscripten_vr_display_presenting(handle)) {
        emscripten_vr_cancel_display_render_loop(handle);
        EM_ASM($('#vr-exit').hide());
        EM_ASM($('#vr-enter').show());
        EM_ASM(resizeCanvas());
    }else{
        VRFrameData data;
        emscripten_vr_get_frame_data(handle, &data);
        vr_update_pos(data);
        gui.drawVR(data.leftProjectionMatrix, data.leftViewMatrix,
                   data.rightProjectionMatrix, data.rightViewMatrix,
                   vrPos, vrWidth, vrHeight);
        emscripten_vr_submit_frame(handle);
    }
}

EM_BOOL vr_stop_presenting(int, const EmscriptenMouseEvent*, void*){
    emscripten_vr_exit_present(handle);
    return 1;
}

EM_BOOL vr_start_presenting(int, const EmscriptenMouseEvent*, void*){
    VRLayerInit layer{NULL, VR_LAYER_DEFAULT_LEFT_BOUNDS, VR_LAYER_DEFAULT_RIGHT_BOUNDS};
    if(!emscripten_vr_request_present(handle, &layer, 1, nullptr, nullptr)){
        printf("Request present failed\n");
        return 0;
    }
    if(!emscripten_vr_set_display_render_loop(handle, vr_loop)){
        printf("Error: failed to set vr-render-loop\n.");
    }
    VREyeParameters leftParams, rightParams;
    emscripten_vr_get_eye_parameters(handle, VREyeLeft, &leftParams);
    emscripten_vr_get_eye_parameters(handle, VREyeRight, &rightParams);
    vrWidth = leftParams.renderWidth = rightParams.renderWidth;
    vrHeight = std::max(leftParams.renderHeight, rightParams.renderHeight);
    emscripten_set_canvas_element_size("#canvas", vrWidth, vrHeight);
    EM_ASM($('#vr-enter').hide());
    EM_ASM($('#vr-exit').show());
    return 1;
}

void tryInitVR(void*){
    printf("Browser running WebVR version %d.%d\n",
                    emscripten_vr_version_major(),
                    emscripten_vr_version_minor());
    if(!emscripten_vr_ready()){
        printf("VR not initialized\n");
        return;
    }else{
        int numDisplays = emscripten_vr_count_displays();
        if(!numDisplays){
            printf("No VR displays connected\n");
            return;
        }
        printf("Number of VR displays: %d\n", numDisplays);
        for(int i=0; i<numDisplays; ++i){
            handle = emscripten_vr_get_display_handle(i);
            printf("Display %d: %s\n", i, emscripten_vr_get_display_name(handle));
            VRDisplayCapabilities caps;
            if(!emscripten_vr_get_display_capabilities(handle, &caps)){
                printf("Error: failed to get display capabilities.\n");
                continue;
            }
            printf("Display Capabilities:\n"
                   "{hasPosition: %d, hasExternalDisplay: %d, canPresent: %d, maxLayers: %u}\n",
                   caps.hasPosition, caps.hasExternalDisplay, caps.canPresent, caps.maxLayers);
            if(caps.hasExternalDisplay && !emscripten_vr_display_connected(handle)){
                printf("Error: has external display, but is not connected.\n");
                continue;
            }
            printf("Succeeded so far, using this display\n");
            vrHasPos = caps.hasPosition;
            EM_ASM($('#vr-enter').show());
            emscripten_set_click_callback("vr-enter", nullptr, 0, vr_start_presenting);
            emscripten_set_click_callback("vr-exit", nullptr, 0, vr_stop_presenting);
            return;
        }
    }
    printf("Could not find a working VR display\n");
}

int main()
{
    // initialize library
    Vipster::readConfig();
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
        EM_ASM(alertWebGL());
        return 0;
    }
    emscripten_vr_init(tryInitVR, nullptr);

    // init GL
    emscripten_webgl_make_context_current(context);
    gui.initGL("# version 300 es\nprecision highp float;\n", "");

    // init examples (something needs to be displayed for the renderer to not fail
    Step* step;
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
    for(auto& fmt: IOPlugins){
        EM_ASM_({addParser($0, $1)}, fmt.first, fmt.second->name.c_str());
        if(fmt.second->writer){
            EM_ASM_({addWriter($0, $1)}, fmt.first, fmt.second->name.c_str());
        }
    }
    EM_ASM(setMol(0));
    emscripten_set_main_loop(main_loop, 0, 1);
    emscripten_webgl_destroy_context(context);
    return 1;
}
