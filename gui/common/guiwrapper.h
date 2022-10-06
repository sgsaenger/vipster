#ifndef GLWRAPPER_H
#define GLWRAPPER_H

#include <string>
#ifdef WEBVIPSTER
#include <GLES3/gl3.h>
#else
#include <QOpenGLExtraFunctions>
#endif
#include <vector>
#include <array>
#include "vipster/molecule.h"
#include "vipster/settings.h"
#include "guimat.h"
#include "guidata.h"
#include "stepdata.h"
#include "seldata.h"
#include "guiglobals.h"

namespace Vipster {

#ifdef WEBVIPSTER
class GuiWrapper{
#else
class GuiWrapper: protected QOpenGLExtraFunctions{
public:
    virtual ~GuiWrapper() = default;
#endif
public:
    GuiWrapper(const Settings &s);
    void initGL();
    void draw(void *context=nullptr);
    void drawSel(void *context=nullptr);
    void drawVR(const float* leftProj, const float* leftView,
                const float* rightProj, const float* rightView,
                const Vec& pos, unsigned long width, unsigned long height);
    void setMainStep(Step* step);
    void setMainSel(Step::selection* sel);
    void updateMainStep();
    void updateMainSelection();
    void resizeViewMat(long w, long h);
    void zoomViewMat(float i);
    void rotateViewMat(float x, float y, float z);
    void translateViewMat(float x, float y, float z);
    enum class alignDir{x,y,z,mx,my,mz};
    void alignViewMat(alignDir d);
    Mat getAxes();
    // cpu-side data
    GUI::PBCVec mult{{1,1,1}};
    Step* curStep{nullptr};
    Step::selection* curSel{nullptr};
    GUI::StepData mainStep{nullptr};
    GUI::SelData selection{nullptr};
    std::vector<std::weak_ptr<GUI::Data>> *stepExtras{nullptr};
    std::vector<std::weak_ptr<GUI::Data>> *vpExtras{nullptr};
    const Settings &settings;
    // gpu-side global data
    GLuint view_ubo;
    // cpu-side global data
    GUI::Mat_16f vMat, oMat, pMat, rMat;
    bool vMatChanged, pMatChanged, rMatChanged;
    // state
    bool drawPerspective{false};
private:
    void drawPre(void *context=nullptr);
    void drawImpl(const Vec& pos, void *context=nullptr);
    void updateViewUBO();
    void updateViewUBOVR(const float* proj, const float* view);
};

}

#endif // GLWRAPPER_H
