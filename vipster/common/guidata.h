#ifndef GUIDATA_H
#define GUIDATA_H

#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#include <string>
#else
#include <QOpenGLFunctions_3_3_Core>
#endif

#include "vec.h"

namespace Vipster{
namespace GUI {

#ifdef __EMSCRIPTEN__
struct GlobalData{
#else
struct GlobalData: protected QOpenGLFunctions_3_3_Core{
#endif
    GlobalData(const std::string& header, const std::string& folder);
    GLuint buffers[3];
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
    virtual ~Data() = default;
    Data(const GlobalData&);
    virtual void drawMol() = 0;
    virtual void drawCell(const std::array<uint8_t,3> &mult) = 0;
    const GlobalData& global;
    virtual void syncToGPU() = 0;
    GLuint loadShader(std::string vert, std::string frag);
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
