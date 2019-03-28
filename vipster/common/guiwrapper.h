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
#include "guiglobals.h"

namespace Vipster {

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
    void drawVR(const float* leftProj, const float* leftView,
                const float* rightProj, const float* rightView,
                const Vec& pos, unsigned long width, unsigned long height);
    void setMainStep(Step* step, bool draw_bonds=true, bool draw_cell=true);
    void setMainSel(Step::selection* sel);
    void updateMainStep(bool draw_bonds=true, bool draw_cell=true);
    void updateMainSelection();
    void addExtraData(GUI::Data* dat);
    void delExtraData(GUI::Data* dat);
public:
    void resizeViewMat(long w, long h);
    void zoomViewMat(float i);
    void rotateViewMat(float x, float y, float z);
    void translateViewMat(float x, float y, float z);
    enum class alignDir{x,y,z,mx,my,mz};
    void alignViewMat(alignDir d);
    Mat getAxes();
    // cpu-side data
    PBCVec mult{{1,1,1}};
    Step* curStep{nullptr};
    Step::selection* curSel{nullptr};
    GUI::GlobalData globals;
    GUI::StepData mainStep{globals, nullptr};
    GUI::SelData selection{globals, settings.selCol.val, nullptr};
    std::set<GUI::Data*> extraData{};
private:
    void drawPre(void);
    void drawImpl(const Vec& pos);
    void updateViewUBO(void);
    void updateViewUBOVR(const float* proj, const float* view);
    // gpu-side global data
    GLuint view_ubo;
    // cpu-side global data
    GUI::Mat vMat, oMat, pMat, rMat;
    bool vMatChanged, pMatChanged, rMatChanged, drawPerspective;
};

}

#endif // GLWRAPPER_H
