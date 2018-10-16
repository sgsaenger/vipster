#ifndef GLWRAPPER_H
#define GLWRAPPER_H

#include <string>
#ifdef __EMSCRIPTEN__
#include <GLES3/gl3.h>
#else
#include <QOpenGLFunctions_3_3_Core>
#endif
#include <vector>
#include <array>
#include "molecule.h"
#include "guimat.h"
#include "guidata.h"
#include "stepdata.h"
#include "seldata.h"

/*
 *
 *
 *
 *  CAREFUL with initialization order of data-classes and guiwrapper!
 *
 *
 */

namespace Vipster {

namespace Enums{
enum GuiChange: uint8_t{
    atoms=0x1, cell=0x2, fmt=0x4, kpoints=0x8,
    selection=0x10, settings=0x20,
    extra=0x40};
}
using Enums::GuiChange;
constexpr auto guiStepChanged = GuiChange::atoms | GuiChange::cell | GuiChange::fmt | GuiChange::selection;
constexpr auto guiMolChanged = GuiChange::kpoints;

#ifdef __EMSCRIPTEN__
class GuiWrapper{
#else
class GuiWrapper: protected QOpenGLFunctions_3_3_Core{
public:
    virtual ~GuiWrapper() = default;
#endif
public:
    void initGL(const std::string& header, const std::string& folder);
    void draw(void);
    void drawSel(void);
    void setMainStep(Step* step, bool draw_bonds=true);
    void setMainSel(Step::selection* sel);
    void updateMainStep(bool draw_bonds=true);
    void updateMainSelection();
    void addExtraData(GUI::Data* dat);
    void delExtraData(GUI::Data* dat);
public:
    void resizeViewMat(int w, int h);
    void zoomViewMat(int i);
    void rotateViewMat(float x, float y, float z);
    void translateViewMat(float x, float y, float z);
    enum class alignDir{x,y,z,mx,my,mz};
    void alignViewMat(alignDir d);
    Mat getAxes();
    // cpu-side data
    std::array<uint8_t,3> mult{{1,1,1}};
    Step* curStep{nullptr};
    Step::selection* curSel{nullptr};
    GUI::GlobalData globals;
    GUI::StepData mainStep{globals, nullptr};
    GUI::SelData selection{globals, nullptr};
    std::set<GUI::Data*> extraData{};
private:
    void updateViewUBO(void);
    // gpu-side global data
    GLuint view_ubo;
    // cpu-side global data
    GUI::Mat vMat, oMat, pMat, rMat;
    bool vMatChanged, pMatChanged, rMatChanged, drawPerspective;
};

}

#endif // GLWRAPPER_H
