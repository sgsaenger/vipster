#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>

#include <iostream>

#include <guiwrapper.h>
#include <molecule.h>
#include <iowrapper.h>
#include <atom_model.h>

namespace em = emscripten;
using namespace Vipster;

static GuiWrapper gui;

template<typename T>
em::class_<std::array<T, 3>> register_array(const char* name) {
    typedef std::array<T, 3> ArrType;
    return em::class_<std::array<T, 3>>(name)
        .template constructor<>()
        .function("size", &ArrType::size)
        .function("get", &em::internal::VectorAccess<ArrType>::get)
        .function("set", &em::internal::VectorAccess<ArrType>::set)
    ;
}

EMSCRIPTEN_BINDINGS(vipster) {
    em::enum_<IOFmt>("IOFmt")
            .value("XYZ", IOFmt::XYZ)
            .value("PWI", IOFmt::PWI)
            ;
}

extern "C" {
EMSCRIPTEN_KEEPALIVE
void emReadFile(char* fn, IOFmt fmt){
    auto d = readFile(fn, fmt);
    gui.molecules.push_back(d.mol);
    gui.updateBuffers(&gui.molecules.back().getStep(0), true);
}

EMSCRIPTEN_KEEPALIVE
void emSetMult(int x, int y, int z){
    gui.mult = {{x,y,z}};
}

}

EM_BOOL mouse_event(int eventType, const EmscriptenMouseEvent* mouseEvent, void* userData)
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
        }
        break;
    }
    localX = mouseEvent->canvasX;
    localY = mouseEvent->canvasY;
    return 1;
}

EM_BOOL touch_event(int eventType, const EmscriptenTouchEvent* touchEvent, void* userData)
{
    constexpr long translateDelta = 10, scaleDelta = 10;
    enum class TwoTouch { None, Scale, Translate};
    static TwoTouch ttMode = TwoTouch::None;
    static long local1X, local2X, local1Y, local2Y, distance, transX, transY;
    long tmp=0, tmp2=0;
    switch (eventType) {
    case EMSCRIPTEN_EVENT_TOUCHSTART:
        switch(touchEvent->numTouches){
        case 2:
            local2X = touchEvent->touches[1].canvasX;
            local2Y = touchEvent->touches[1].canvasY;
            distance = std::sqrt(std::pow(local2X-touchEvent->touches[0].canvasX,2)+
                                 std::pow(local2Y-touchEvent->touches[0].canvasY,2));
            transX = (local2X + touchEvent->touches[0].canvasX)/2;
            transY = (local2Y + touchEvent->touches[0].canvasY)/2;
        case 1:
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
            gui.rotateViewMat(touchEvent->touches[0].canvasX-local1X,
                    touchEvent->touches[0].canvasY-local1Y,0);
            local1X = touchEvent->touches[0].canvasX;
            local1Y = touchEvent->touches[0].canvasY;
            break;
        case 2:
            if(ttMode == TwoTouch::None){
                tmp = std::sqrt(std::pow(touchEvent->touches[1].canvasX-
                                         touchEvent->touches[0].canvasX,2)+
                                std::pow(touchEvent->touches[1].canvasY-
                                         touchEvent->touches[0].canvasY,2));
                if (std::abs(tmp-distance)>scaleDelta){
                    ttMode = TwoTouch::Scale;
                }else if(local1X-touchEvent->touches[0].canvasX > translateDelta||
                         local1Y-touchEvent->touches[0].canvasY > translateDelta||
                         local2X-touchEvent->touches[1].canvasX > translateDelta||
                         local2Y-touchEvent->touches[1].canvasY > translateDelta){
                    ttMode = TwoTouch::Translate;
                }
            }
            if(ttMode == TwoTouch::None) break;
            if(ttMode == TwoTouch::Scale){
                if(!tmp)
                    tmp = std::sqrt(std::pow(touchEvent->touches[1].canvasX-
                                             touchEvent->touches[0].canvasX,2)+
                                    std::pow(touchEvent->touches[1].canvasY-
                                             touchEvent->touches[0].canvasY,2));
                gui.zoomViewMat(tmp-distance);
                distance = tmp;
            }else if(ttMode == TwoTouch::Translate){
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
        ttMode = TwoTouch::None;
        break;
    }
    return 1;
}

EM_BOOL wheel_event(int eventType, const EmscriptenWheelEvent* wheelEvent, void* userData)
{
    gui.zoomViewMat(-wheelEvent->deltaY);
    return 1;
}

EM_BOOL key_event(int eventType, const EmscriptenKeyboardEvent* keyEvent, void* userData)
{
    switch(keyEvent->keyCode){
    case 37:
        gui.rotateViewMat(-10, 0, 0);
        break;
    case 39:
        gui.rotateViewMat(10, 0, 0);
        break;
    case 38:
        gui.rotateViewMat(0, -10, 0);
        break;
    case 40:
        gui.rotateViewMat(0, 10, 0);
        break;
    }
    return 1;
}

void one_iter(){
    int width, height, fullscreen;
    static int localWidth, localHeight;
    // handle resize
    emscripten_get_canvas_size(&width, &height, &fullscreen);
    if( width != localWidth || height != localHeight){
        gui.resizeViewMat(width, height);
        localWidth = width;
        localHeight = height;
    }
    // sync data and draw
    gui.updateVBOs();
    gui.updateUBO();
    gui.draw();
}

int main()
{
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
    gui.initUBO();
    gui.initViewMat();

    gui.molecules.emplace_back("test");
    Step* step = &gui.molecules[0].getStep(0);
    step->newAtom();
    step->newAtom({"O",{{1,0,0}}});
    step->newAtom({"F",{{0,1,0}}});
    gui.updateBuffers(step, true);

    emscripten_set_mousedown_callback("#canvas", nullptr, 1, mouse_event);
    emscripten_set_mouseup_callback(0, nullptr, 1, mouse_event);
    emscripten_set_mousemove_callback(0, nullptr, 1, mouse_event);
    emscripten_set_wheel_callback("#canvas", nullptr, 1, wheel_event);
    emscripten_set_keypress_callback(0, nullptr, 1, key_event);
    emscripten_set_touchstart_callback("#canvas", nullptr, 1, touch_event);
    emscripten_set_touchmove_callback("#canvas", nullptr, 1, touch_event);
    emscripten_set_touchend_callback("#canvas", nullptr, 1, touch_event);
    emscripten_set_main_loop(one_iter, 0, 1);
    gui.deleteGLObjects();
    return 1;
}
