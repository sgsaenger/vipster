#include <iostream>
#include <sys/types.h>
#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>
#include <GLES3/gl3.h>
#include <molecule.h>
#include <atom_model.h>

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
    register_array<int>("intarr");
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

    GLuint atom_prog = glCreateProgram();
    GLuint atom_vs = glCreateShader(GL_VERTEX_SHADER);
    GLuint atom_fs = glCreateShader(GL_FRAGMENT_SHADER);
    int gl_ok = 0;

    const char *atomVert =
        "#version 300 es\n"
        "in vec3 vertex_modelspace;\n"
        "in vec3 position_modelspace;\n"
        "in float scale_modelspace;\n"
        "in vec4 color_input;\n"
        "uniform mat4 vpMatrix;\n"
        "uniform mat4 rMatrix;\n"
        "uniform float atom_fac;\n"
        "uniform vec3 offset;\n"
        "out vec3 normals_cameraspace;\n"
        "out vec3 EyeDirection_cameraspace;\n"
        "out vec3 LightDirection_cameraspace;\n"
        "out vec4 MaterialDiffuseColor;\n"
        "void main(void)\n"
        "{\n"
        "   gl_Position = vpMatrix * vec4(vertex_modelspace * scale_modelspace *"
        " atom_fac + position_modelspace + offset, 1);\n"
        "   vec3 vertex_cameraspace = (rMatrix * vec4(vertex_modelspace, 1)).xyz;\n"
        "   EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;\n"
        "   LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;\n"
        "   normals_cameraspace = (rMatrix * vec4(vertex_modelspace, 0)).xyz;\n"
        "   MaterialDiffuseColor = color_input;\n"
        "}\n";

    const char *atomFrag = 
        "# version 300 es\n"
        "precision highp float;\n"
        "out vec4 fragColor;\n"
        "in vec3 normals_cameraspace;\n"
        "in vec3 EyeDirection_cameraspace;\n"
        "in vec3 LightDirection_cameraspace;\n"
        "in vec4 MaterialDiffuseColor;\n"
        "void main(void)\n"
        "{\n"
        "   vec3 LightColor = vec3(1,1,1);\n"
        "   float LightPower = 50.f;\n"
        "   float distance = 10.;\n"
        "   vec3 MaterialAmbientColor = vec3(0.5,0.5,0.5) * MaterialDiffuseColor.xyz;\n"
        "   vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);\n"
        "   vec3 n = normalize(normals_cameraspace);\n"
        "   vec3 l = normalize(LightDirection_cameraspace);\n"
        "   float cosTheta = clamp(dot(n,l),0.,1.);\n"
        "   vec3 E = normalize(EyeDirection_cameraspace);\n"
        "   vec3 R = reflect(-l,n);\n"
        "   float cosAlpha = clamp(dot(E,R),0.,1.);\n"
        "   vec3 fragTemp = MaterialAmbientColor +\n"
        "MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,10.)/(distance*distance) + \n"
        "MaterialDiffuseColor.xyz * LightColor * LightPower * cosTheta/(distance*distance);\n"
        "   fragColor = vec4(fragTemp,MaterialDiffuseColor.a);\n"
        "}\n";

    glShaderSource(atom_vs, 1, &atomVert, NULL);
    glCompileShader(atom_vs);
    glGetShaderiv(atom_vs, GL_COMPILE_STATUS, &gl_ok);
    if(!gl_ok){
        GLint infoLen = 0;
        glGetShaderiv(atom_vs, GL_INFO_LOG_LENGTH, &infoLen);
        if(infoLen > 1)
        {
            char* infoLog = (char*) malloc(sizeof(char)*infoLen+1);
            glGetShaderInfoLog(atom_vs,infoLen,NULL,infoLog);
            printf("Error compiling vertshader:\n%s\n", infoLog);
            free(infoLog);
        }
    }else{
        printf("vert-shader sucess\n");
    }
    glShaderSource(atom_fs, 1, &atomFrag, NULL);
    glCompileShader(atom_fs);
    glGetShaderiv(atom_fs, GL_COMPILE_STATUS, &gl_ok);
    if(!gl_ok){
        printf("frag-shader failure\n");
        GLint infoLen = 0;
        glGetShaderiv(atom_fs, GL_INFO_LOG_LENGTH, &infoLen);
        if(infoLen > 1)
        {
            char* infoLog = (char*) malloc(sizeof(char)*infoLen+1);
            glGetShaderInfoLog(atom_fs,infoLen,NULL,infoLog);
            printf("Error compiling fragshader:\n%s\n", infoLog);
            free(infoLog);
        }
    }else{
        printf("frag-shader sucess\n");
    }

    glAttachShader(atom_prog,atom_vs);
    glAttachShader(atom_prog,atom_fs);
    glLinkProgram(atom_prog);
    glGetProgramiv(atom_prog, GL_LINK_STATUS, &gl_ok);
    if(!gl_ok){
        printf("linking error\n");
    }
    glUseProgram(atom_prog);

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
    GLint vpMatLoc = glGetUniformLocation(atom_prog, "vpMatrix");
    glUniformMatrix4fv(vpMatLoc, 1, false, vpMat);
    GLint rMatLoc = glGetUniformLocation(atom_prog, "rMatrix");
    glUniformMatrix4fv(rMatLoc, 1, false, rMat);
    GLint atfacLoc = glGetUniformLocation(atom_prog, "atom_fac");
    glUniform1f(atfacLoc, (GLfloat)0.5);
    GLint offsetLoc = glGetUniformLocation(atom_prog, "offset");
    glUniform3fv(offsetLoc, 1, offset);

    Molecule mol{"test"};
    Step& st = mol.getStep(0);
    st.newAtom();
    st.newAtom({"O",{1,0,0}});
    st.newAtom({"F",{0,1,0}});
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

    return 1;
}
