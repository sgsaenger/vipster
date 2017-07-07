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
    gui.curStep = &gui.molecules.back().getStep(0);
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
                gui.rMat = {{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}};
                gui.rMatChanged = true;
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
            guiMatRot(gui.rMat, mouseEvent->canvasX-localX, 0, 1, 0);
            guiMatRot(gui.rMat, mouseEvent->canvasY-localY, 1, 0, 0);
            gui.rMatChanged = true;
            break;
        case OpMode::Translation:
            guiMatTranslate(gui.vMat, (mouseEvent->canvasX-localX)/10.,
                           -(mouseEvent->canvasY-localY)/10.,0);
            gui.vMatChanged = true;
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
            gui.rMat = {{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}};
            gui.rMatChanged = true;
            break;
        }
        break;
    case EMSCRIPTEN_EVENT_TOUCHMOVE:
        switch(touchEvent->numTouches){
        case 1:
            guiMatRot(gui.rMat, touchEvent->touches[0].canvasX-local1X, 0, 1, 0);
            guiMatRot(gui.rMat, touchEvent->touches[0].canvasY-local1Y, 1, 0, 0);
            gui.rMatChanged = true;
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
                guiMatScale(gui.vMat, (tmp-distance)<0?1.1:0.9);
                gui.vMatChanged = true;
                distance = tmp;
            }else if(ttMode == TwoTouch::Translate){
                tmp = (touchEvent->touches[0].canvasX + touchEvent->touches[1].canvasX)/2;
                tmp2 = (touchEvent->touches[0].canvasY + touchEvent->touches[1].canvasY)/2;
                guiMatTranslate(gui.vMat, (tmp-transX)/10., (tmp-transY)/10., 0);
                gui.vMatChanged = true;
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
    guiMatScale(gui.vMat, wheelEvent->deltaY<0?1.1:0.9);
    gui.vMatChanged = true;
    return 1;
}

EM_BOOL key_event(int eventType, const EmscriptenKeyboardEvent* keyEvent, void* userData)
{
    switch(keyEvent->keyCode){
        case 37:
            guiMatRot(gui.rMat, -10, 0, 1, 0);
            gui.rMatChanged = true;
        break;
        case 39:
            guiMatRot(gui.rMat, 10, 0, 1, 0);
            gui.rMatChanged = true;
        break;
        case 38:
            guiMatRot(gui.rMat, -10, 1, 0, 0);
            gui.rMatChanged = true;
        break;
        case 40:
            guiMatRot(gui.rMat, 10, 1, 0, 0);
            gui.rMatChanged = true;
        break;
    }
    return 1;
}

void one_iter(){
    int width, height, fullscreen;
    static const Step* step{nullptr};
    static int localWidth, localHeight;
    // update Step-data
    if(step != gui.curStep){
        step = gui.curStep;
        gui.atom_buffer.clear();
        gui.atom_buffer.reserve(gui.curStep->getNat());
        for(const Atom& at:gui.curStep->getAtoms()){
            PseEntry &pse = (*gui.curStep->pse)[at.name];
            gui.atom_buffer.push_back({{at.coord[0],at.coord[1],at.coord[2],pse.covr,
                                    pse.col[0],pse.col[1],pse.col[2],pse.col[3]}});
        }
        glBindBuffer(GL_ARRAY_BUFFER, gui.atom_vbo);
        glBufferData(GL_ARRAY_BUFFER, gui.atom_buffer.size()*8*sizeof(float),
                     (void*)gui.atom_buffer.data(), GL_STREAM_DRAW);
    }
    // handle resize
    emscripten_get_canvas_size(&width, &height, &fullscreen);
    if( width != localWidth || height != localHeight){
        height==0?height=1:0;
        glViewport(0,0,width,height);
        float aspect = (float)width/height;
        gui.pMat = guiMatMkOrtho(-10*aspect,10*aspect,-10, 10, 0, 1000);
        gui.pMatChanged = true;
        localWidth = width;
        localHeight = height;
    }
    // update mvp-matrices if needed
    if(gui.rMatChanged){
        glUniformMatrix4fv(glGetUniformLocation(gui.atom_program, "rMatrix"),
                           1, true, gui.rMat.data());
        glUniformMatrix4fv(glGetUniformLocation(gui.atom_program, "vpMatrix"),
                           1, true, (gui.pMat*gui.vMat*gui.rMat).data());
        gui.pMatChanged = gui.vMatChanged = gui.rMatChanged = false;
    }else if(gui.pMatChanged || gui.vMatChanged){
        glUniformMatrix4fv(glGetUniformLocation(gui.atom_program, "vpMatrix"),
                           1, true, (gui.pMat*gui.vMat*gui.rMat).data());
        gui.pMatChanged = gui.vMatChanged = false;
    }
    // draw stuff
    glBindVertexArray(gui.atom_vao);
    glDrawArraysInstanced(GL_TRIANGLES,0,atom_model_npoly,gui.atom_buffer.size());
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

    gui.loadShader(gui.atom_program, "# version 300 es\nprecision highp float;\n",
                   readShader("/atom.vert"),
                   readShader("/atom.frag"));
    glUseProgram(gui.atom_program);

    gui.initAtomVAO();

    gui.vMat = guiMatMkLookAt({{0,0,10}},{{0,0,0}},{{0,1,0}});
    gui.rMat = {{1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1}};
    gui.pMatChanged = gui.vMatChanged = gui.rMatChanged = true;
    Vec offset = {{0,0,0}};
    GLuint atfacLoc = glGetUniformLocation(gui.atom_program, "atom_fac");
    glUniform1f(atfacLoc, (GLfloat)0.5);
    GLuint offsetLoc = glGetUniformLocation(gui.atom_program, "offset");
    glUniform3fv(offsetLoc, 1, offset.data());

    gui.molecules.emplace_back("test");
    gui.curStep = &gui.molecules[0].getStep(0);
    Step* step = const_cast<Step*>(gui.curStep);
    step->newAtom();
    step->newAtom({"O",{{1,0,0}}});
    step->newAtom({"F",{{0,1,0}}});


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
