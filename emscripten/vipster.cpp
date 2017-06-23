#include <sys/types.h>
#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>
#include <GLES3/gl3.h>
#include <molecule.h>
#include <atom_model.h>
#include <glwrapper.h>

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

void one_iter(){
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

    GLuint atomProg = loadShader("# version 300 es\nprecision highp float;\n",
                                 "/atom.vert",
                                 "/atom.frag");
    glUseProgram(atomProg);

    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

    GLuint sphere_vbo;
    glGenBuffers(1,&sphere_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    glBufferData(GL_ARRAY_BUFFER, atom_model_npoly*3*sizeof(float), (void*)&atom_model, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, NULL);

    float vpMat[16] = {1,0,0,0,
                       0,1,0,0,
                       0,0,1,0,
                       0,0,0,1};
    float rMat[16] = {1,0,0,0,
                      0,1,0,0,
                      0,0,1,0,
                      0,0,0,1};
    float offset[3] = {0,0,0};
    GLint vpMatLoc = glGetUniformLocation(atomProg, "vpMatrix");
    glUniformMatrix4fv(vpMatLoc, 1, false, vpMat);
    GLint rMatLoc = glGetUniformLocation(atomProg, "rMatrix");
    glUniformMatrix4fv(rMatLoc, 1, false, rMat);
    GLint atfacLoc = glGetUniformLocation(atomProg, "atom_fac");
    glUniform1f(atfacLoc, (GLfloat)0.5);
    GLint offsetLoc = glGetUniformLocation(atomProg, "offset");
    glUniform3fv(offsetLoc, 1, offset);

    Molecule mol{"test"};
    Step& st = mol.getStep(0);
    st.newAtom();
    st.newAtom({"O",{{1,0,0}}});
    st.newAtom({"F",{{0,1,0}}});
    std::vector<std::array<float,8>> atom_buffer;
    for(const Atom& at:st.getAtoms()){
        PseEntry &pse = (*st.pse)[at.name];
        atom_buffer.push_back({{at.coord[0],at.coord[1],at.coord[2],pse.covr,
                                pse.col[0],pse.col[1],pse.col[2],pse.col[3]}});
    }
    GLuint atom_vbo;
    glGenBuffers(1,&atom_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
    glBufferData(GL_ARRAY_BUFFER, atom_buffer.size()*8*sizeof(float), (void*)atom_buffer.data(), GL_STATIC_DRAW);

    //DRAW!
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,8*sizeof(float),0);
    glVertexAttribDivisor(1,1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,8*sizeof(float),(const GLvoid*)(3*sizeof(float)));
    glVertexAttribDivisor(2,1);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,8*sizeof(float),(const GLvoid*)(4*sizeof(float)));
    glVertexAttribDivisor(3,1);

    glDrawArraysInstanced(GL_TRIANGLES,0,atom_model_npoly,atom_buffer.size());

    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDisableVertexAttribArray(2);
    glDisableVertexAttribArray(3);
    glVertexAttribDivisor(1,0);
    glVertexAttribDivisor(2,0);
    glVertexAttribDivisor(3,0);
    glBindBuffer(GL_ARRAY_BUFFER, NULL);

    emscripten_set_main_loop(one_iter, 0, 1);
    return 1;
}
