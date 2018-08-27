#ifndef GUIDATA_H
#define GUIDATA_H

#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
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
    QOpenGLContext* context;
#endif
    GlobalData(const std::string& header, const std::string& folder);
    GLuint atom_program, bond_program, cell_program, sel_program;
    GLuint buffers[3];
    GLuint& sphere_vbo, &cylinder_vbo;
    GLuint& cell_ibo;
private:
    void loadShader(GLuint &program, const std::string &header,
                    std::string vertShaderStr, std::string fragShaderStr);
};

#ifdef __EMSCRIPTEN__
class Data
#else
class Data: protected QOpenGLFunctions_3_3_Core
#endif
{
public:
    virtual ~Data() = default;
    Data(GlobalData&);
    virtual void drawMol() = 0;
    virtual void drawCell(const std::array<uint8_t,3> &mult) = 0;
    GlobalData& global;
    virtual void syncToGPU() = 0;
};

}
}

#endif // GUIDATA_H
