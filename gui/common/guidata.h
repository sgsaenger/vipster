#ifndef GUIDATA_H
#define GUIDATA_H

#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#include <string>
#include <vector>
#include <map>
#else
#include <QOpenGLExtraFunctions>
#endif

#include "vipster/vec.h"
#include "guiglobals.h"

namespace Vipster{
namespace GUI {

/* Global OpenGL Data
 *
 * Contains strings to find/compile shaders correctly (platform/runtime dependent)
 * indices for shared OpenGL objects like basic mesh data
 */
#ifdef __EMSCRIPTEN__
struct GlobalData{
#else
struct GlobalData: protected QOpenGLExtraFunctions{
#endif
    GlobalData();
    void initGL();
    GLuint sphere_vbo, cylinder_vbo;
    GLuint cell_ibo;
    std::string header, folder;
    bool initialized{false};
};

/* Base for concrete OpenGL Wrappers
 *
 * provides mechanism to load shaders
 * public interface for syncing data between cpu/gpu
 */
#ifdef __EMSCRIPTEN__
class Data
#else
class Data: protected QOpenGLExtraFunctions
#endif
{
public:
    Data(const GlobalData&);
    virtual ~Data() = default;
    Data(Data&&);
    Data(const Data&) = delete;
    Data& operator=(const Data&) = delete;
    Data& operator=(Data&&) = delete;
    const GlobalData& global;

    virtual void draw(const Vec &off, const PBCVec &mult,
                      const Mat &cv, bool drawCell, void *context) = 0;
    GLuint loadShader(const std::string &vert, const std::string &frag);
    // TODO: split this in shader/vbo/vao management to reduce context sensitivity
    void syncToGPU(void *context);

protected:
    virtual void updateGL() = 0;
    virtual void initGL(void *context) = 0;
    std::map<void*, bool> initialized{};
    bool updated{true};
#ifndef __EMSCRIPTEN__
    bool wrap_initialized{false};
#endif
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
