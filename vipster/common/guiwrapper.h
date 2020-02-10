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
    GuiWrapper(GUI::GlobalData &g, const Settings &s);
    void initGL(const std::string& header, const std::string& folder);
    void draw(void *context=nullptr);
    void drawSel(void *context=nullptr);
    void drawVR(const float* leftProj, const float* leftView,
                const float* rightProj, const float* rightView,
                const Vec& pos, unsigned long width, unsigned long height);
    void setMainStep(Step* step);
    void setMainSel(Step::selection* sel);
    void updateMainStep();
    void updateMainSelection();
public:
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
    GUI::GlobalData& globals;
    GUI::StepData mainStep{globals, nullptr};
    GUI::SelData selection{globals, nullptr};
    std::vector<std::shared_ptr<GUI::Data>> *stepExtras{nullptr};
    std::vector<std::shared_ptr<GUI::Data>> *vpExtras{nullptr};
private:
    const Settings &settings;
    void drawPre(void *context=nullptr);
    void drawImpl(const Vec& pos, void *context=nullptr);
    void updateViewUBO();
    void updateViewUBOVR(const float* proj, const float* view);
    // gpu-side global data
    GLuint view_ubo;
    // cpu-side global data
    /* rotation matrix
     * converts from model space to world space
     * rotation only as shifting is done by the data itself
     */
    GUI::Mat_16f rMat;
    bool rMatChanged;
    /* view matrix
     * converts world space to view space
     */
    GUI::Mat_16f vMat;
    bool vMatChanged;
    /* projection matrices
     * converts view space to clip space
     * precomputed for orthogonal and perspective projection
     */
    GUI::Mat_16f oMat, pMat;
    bool pMatChanged;
    // state
    bool drawPerspective{false};
};

}

#endif // GLWRAPPER_H
