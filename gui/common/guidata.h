#ifndef GUIDATA_H
#define GUIDATA_H

#ifdef WEBVIPSTER
#include <GLES3/gl3.h>
#else
#include <QOpenGLExtraFunctions>
#endif
#include <string>
#include <vector>
#include <map>

#include "vipster/vec.h"
#include "guiglobals.h"

namespace Vipster{
namespace GUI {

/* Base for concrete OpenGL Wrappers
 *
 * provides mechanism to load shaders
 * public interface for syncing data between cpu/gpu
 * manages global OpenGL state for a context
 */
#ifdef WEBVIPSTER
class Data
#else
class Data: protected QOpenGLExtraFunctions
#endif
{
public:
    Data() = default;
    virtual ~Data() = default;
    Data(Data&&);
    Data(const Data&) = delete;
    Data& operator=(const Data&) = delete;
    Data& operator=(Data&&) = delete;

    virtual void draw(const Vec &off, const PBCVec &mult,
                      const Mat &cv, bool drawCell, void *context) = 0;
    void syncToGPU(void *context);

protected:
    struct GlobalContext{
        bool initialized{false};
        GLuint sphere_vbo{0}, cylinder_vbo{0};
        GLuint cell_ibo{0};
        std::string header{}, folder{};
    };
    static std::map<void*, GlobalContext> global_map;
    struct InstanceContext{
        bool initialized{false};
        bool synchronized{false};
    };
    std::map<void*, InstanceContext>  instance_map;
#ifndef WEBVIPSTER
    bool wrap_initialized{false};
#endif
    void initGlobal(void *context);
    virtual void updateGL(void *context) = 0;
    virtual void initGL(void *context) = 0;
    GLuint loadShader(const GlobalContext& globals, const std::string &vert, const std::string &frag);
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
