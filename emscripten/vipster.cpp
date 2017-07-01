#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>
#include <GLES3/gl3.h>
#include <molecule.h>
#include <atom_model.h>
#include <glwrapper.h>
#include <iostream>

namespace em = emscripten;
using namespace Vipster;

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
}

extern "C" {
EMSCRIPTEN_KEEPALIVE
void resizeGL(int w, int h)
{
//    h==0?h=1:0;
    glViewport(0,0,w,h);
//    float aspect = float(w)/h;
//    pMatrix.setToIdentity();
//    //pMatrix.perspective(60.0,aspect,0.001,1000);
//    pMatrix.ortho(-10*aspect,10*aspect,-10,10,0,1000);
}
}

EM_BOOL mouse_event(int eventType, const EmscriptenMouseEvent* mouseEvent, void* userData)
{
    static bool handling_in_progress=false;
    static long localX, localY;
    auto* gl_p = (GLWrapper*)userData;
    switch (eventType) {
    case EMSCRIPTEN_EVENT_MOUSEDOWN:
        if(!handling_in_progress){
            handling_in_progress = true;
            localX = mouseEvent->canvasX;
            localY = mouseEvent->canvasY;
        }
        break;
    case EMSCRIPTEN_EVENT_MOUSEUP:
        if(handling_in_progress){
            handling_in_progress = false;
        }
        break;
    case EMSCRIPTEN_EVENT_MOUSEMOVE:
        if(handling_in_progress){
            glMatRot(gl_p->rMat, mouseEvent->canvasX-localX, 0, 1, 0);
            glMatRot(gl_p->rMat, mouseEvent->canvasY-localY, 1, 0, 0);
            gl_p->rMatChanged = true;
            localX = mouseEvent->canvasX;
            localY = mouseEvent->canvasY;
        }
        break;
    }
    return 1;
}

EM_BOOL wheel_event(int eventType, const EmscriptenWheelEvent* wheelEvent, void* userData)
{
    auto* gl_p = (GLWrapper*)userData;
    glMatScale(gl_p->vMat, wheelEvent->deltaY<0?1.1:0.9);
    gl_p->vMatChanged = true;
    return 1;
}

EM_BOOL key_event(int eventType, const EmscriptenKeyboardEvent* keyEvent, void* userData)
{
    auto* gl_p = (GLWrapper*)userData;
    switch(keyEvent->keyCode){
        case 37:
            glMatRot(gl_p->rMat, -10, 0, 1, 0);
            gl_p->rMatChanged = true;
        break;
        case 39:
            glMatRot(gl_p->rMat, 10, 0, 1, 0);
            gl_p->rMatChanged = true;
        break;
        case 38:
            glMatRot(gl_p->rMat, -10, 1, 0, 0);
            gl_p->rMatChanged = true;
        break;
        case 40:
            glMatRot(gl_p->rMat, 10, 1, 0, 0);
            gl_p->rMatChanged = true;
        break;
    }
    return 1;
}

std::ostream& operator<<(std::ostream &s, const glMat &d)
{
    s << "glMat:\n"
      << d[0] << ' ' << d[1] << ' ' << d[2] << ' ' << d[3] << '\n'
      << d[4] << ' ' << d[5] << ' ' << d[6] << ' ' << d[7] << '\n'
      << d[8] << ' ' << d[9] << ' ' << d[10] << ' ' << d[11] << '\n'
      << d[12] << ' ' << d[13] << ' ' << d[14] << ' ' << d[15] << '\n';
    return s;
}

void one_iter(void *gl_v){
    auto* gl_p = (GLWrapper*)gl_v;
    if(gl_p->rMatChanged){
        glUniformMatrix4fv(glGetUniformLocation(gl_p->atom_program, "rMatrix"),
                           1, true, gl_p->rMat.data());
        glUniformMatrix4fv(glGetUniformLocation(gl_p->atom_program, "vpMatrix"),
                           1, true, (gl_p->pMat*gl_p->vMat*gl_p->rMat).data());
        gl_p->pMatChanged = false;
        gl_p->vMatChanged = false;
        gl_p->rMatChanged = false;
    }else if(gl_p->pMatChanged || gl_p->vMatChanged){
        glUniformMatrix4fv(glGetUniformLocation(gl_p->atom_program, "vpMatrix"),
                           1, true, (gl_p->pMat*gl_p->vMat*gl_p->rMat).data());
        gl_p->pMatChanged = false;
        gl_p->vMatChanged = false;
    }
    glBindVertexArray(gl_p->atom_vao);
    glDrawArraysInstanced(GL_TRIANGLES,0,atom_model_npoly,gl_p->atom_buffer.size());
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
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

    GLWrapper gl{};

    gl.atom_program = loadShader("# version 300 es\nprecision highp float;\n",
                                 "/atom.vert",
                                 "/atom.frag");
    glUseProgram(gl.atom_program);

    float aspect = 16/9;
    gl.pMat = glMatOrtho(-10*aspect,10*aspect,-10,10,0,1000);
//    gl.pMat = {{1,0,0,0,
//                0,1,0,0,
//                0,0,1,0,
//                0,0,0,1}};
    gl.vMat = glMatLookAt({{0,0,10}},{{0,0,0}},{{0,1,0}});
//    gl.vMat = {{1,0,0,0,
//                0,1,0,0,
//                0,0,1,0,
//                0,0,0,1}};
    gl.rMat = {{1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1}};
    gl.offset = {{0,0,0}};
//    GLint vpMatLoc = glGetUniformLocation(gl.atom_program, "vpMatrix");
//    glUniformMatrix4fv(vpMatLoc, 1, false, (gl.pMat*gl.vMat*gl.rMat).data());
//    GLint rMatLoc = glGetUniformLocation(gl.atom_program, "rMatrix");
//    glUniformMatrix4fv(rMatLoc, 1, false, gl.rMat.data());
    GLint atfacLoc = glGetUniformLocation(gl.atom_program, "atom_fac");
    glUniform1f(atfacLoc, (GLfloat)0.5);
    GLint offsetLoc = glGetUniformLocation(gl.atom_program, "offset");
    glUniform3fv(offsetLoc, 1, gl.offset.data());
    gl.pMatChanged = true;
    gl.vMatChanged = true;
    gl.rMatChanged = true;
    gl.offsetChanged = false;

    Molecule mol{"test"};
    Step& st = mol.getStep(0);
    st.newAtom();
    st.newAtom({"O",{{1,0,0}}});
    st.newAtom({"F",{{0,1,0}}});
    for(const Atom& at:st.getAtoms()){
        PseEntry &pse = (*st.pse)[at.name];
        gl.atom_buffer.push_back({{at.coord[0],at.coord[1],at.coord[2],pse.covr,
                                pse.col[0],pse.col[1],pse.col[2],pse.col[3]}});
    }
    glBindBuffer(GL_ARRAY_BUFFER, gl.atom_vbo);
    glBufferData(GL_ARRAY_BUFFER, gl.atom_buffer.size()*8*sizeof(float), (void*)gl.atom_buffer.data(), GL_STATIC_DRAW);


    emscripten_set_mousedown_callback("#canvas", &gl, 1, mouse_event);
    emscripten_set_mouseup_callback(0, &gl, 1, mouse_event);
    emscripten_set_mousemove_callback(0, &gl, 1, mouse_event);
    emscripten_set_wheel_callback("#canvas", &gl, 1, wheel_event);
    emscripten_set_keypress_callback(0, &gl, 1, key_event);
    emscripten_set_main_loop_arg(one_iter, &gl, 0, 1);
    return 1;
}
