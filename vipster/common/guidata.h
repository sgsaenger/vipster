#ifndef GUIDATA_H
#define GUIDATA_H

#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#include <string>
#include <vector>
#else
#include <QOpenGLFunctions_3_3_Core>
#endif

#include "vec.h"
#include "guiglobals.h"

namespace Vipster{
namespace GUI {

#ifdef __EMSCRIPTEN__
struct GlobalData{
#else
struct GlobalData: protected QOpenGLFunctions_3_3_Core{
#endif
    GlobalData();
    void initGL(const std::string& header, const std::string& folder);
    GLuint buffers[3]{0,0,0};
    GLuint& sphere_vbo, &cylinder_vbo;
    GLuint& cell_ibo;
    std::string header, folder;
};

#ifdef __EMSCRIPTEN__
class Data
#else
class Data: protected QOpenGLFunctions_3_3_Core
#endif
{
public:
    const GlobalData& global;

    virtual void drawMol(const Vec &off) = 0;
    virtual void drawCell(const Vec &off, const PBCVec &mult) = 0;
    virtual void updateGL() = 0;
    virtual void initGL() = 0;
    GLuint loadShader(std::string vert, std::string frag);
    void syncToGPU();

    bool updated{true}, initialized{false};
#ifndef __EMSCRIPTEN__
    bool wrap_initialized{false};
#endif

    Data(const GlobalData&);
    virtual ~Data() = default;
    Data(Data&&);
    Data(const Data&) = delete;
    Data& operator=(const Data&) = delete;
    Data& operator=(Data&&) = delete;
};

#define READATTRIB(shader, name) \
    {GLint tmp = glGetAttribLocation(shader.program, #name); \
     if(tmp<0){ \
         throw Vipster::Error("Shader attribute mismatch: "#shader"."#name); \
     }else{ \
         shader.name = static_cast<GLuint>(tmp); \
     }\
    }

#define READUNIFORM(shader, name) \
    {GLint tmp = glGetUniformLocation(shader.program, #name); \
     if(tmp<0){ \
         throw Vipster::Error("Shader uniform mismatch: "#shader"."#name); \
     }else{ \
         shader.name = tmp; \
     }\
    }

}
}

#endif // GUIDATA_H
