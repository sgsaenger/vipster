#include <limits>

#include "guiwrapper.h"

using namespace Vipster;

void GuiWrapper::initGL(const std::string& header, const std::string& folder)
{
#ifndef __EMSCRIPTEN__
    initializeOpenGLFunctions();
#endif
    glClearColor(1,1,1,1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // init globals
    globals.initGL(header, folder);

    // init ViewUBO
    glGenBuffers(1, &view_ubo);
    glBindBuffer(GL_UNIFORM_BUFFER, view_ubo);
    glBufferData(GL_UNIFORM_BUFFER, 2*sizeof(GUI::Mat), nullptr, GL_STATIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, view_ubo);

    // init view matrices
    oMat = guiMatMkOrtho(-15, 15, -10, 10, -100, 1000);
    pMat = guiMatMkPerspective(60.0, 1.5, 0.001f, 1000);
    rMat = {{1,0,0,0,
             0,1,0,0,
             0,0,1,0,
             0,0,0,1}};
    vMat = guiMatMkLookAt({{0,0,10}}, {{0,0,0}}, {{0,1,0}});
    pMatChanged = rMatChanged = vMatChanged = true;
    drawPerspective = settings.perspective.val;
}

void GuiWrapper::draw(void)
{
#ifndef __EMSCRIPTEN__
    if(settings.antialias.val){
        glEnable(GL_MULTISAMPLE);
    }else{
        glDisable(GL_MULTISAMPLE);
    }
#endif
    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(drawPerspective != settings.perspective.val){
        drawPerspective = settings.perspective.val;
        pMatChanged = true;
    }
    if(rMatChanged||pMatChanged||vMatChanged){
        updateViewUBO();
    }
    mainStep.syncToGPU();
    selection.syncToGPU();
    for(const auto& i: extraData){
        i->syncToGPU();
    }
    Vec off = -curStep->getCenter(CdmFmt::Bohr, settings.rotCom.val);
    if(curStep->hasCell()){
        Mat cv = curStep->getCellVec() *
                 curStep->getCellDim(CdmFmt::Bohr);
        off -= (mult[0]-1)*cv[0]/2.;
        off -= (mult[1]-1)*cv[1]/2.;
        off -= (mult[2]-1)*cv[2]/2.;
        mainStep.drawCell(off, mult);
        selection.drawCell(off, mult);
        for(const auto& i: extraData){
            i->drawCell(off, mult);
        }
    }else{
        mainStep.drawMol(off);
        selection.drawMol(off);
        for(const auto& i: extraData){
            i->drawMol(off);
        }
    }
}

void GuiWrapper::drawSel()
{
    mainStep.updateGL();
    mainStep.drawSel(mult);
}

void GuiWrapper::setMainStep(Step* step, bool draw_bonds, bool draw_cell)
{
    curStep = step;
    mainStep.update(step, draw_bonds, draw_cell);
}

void GuiWrapper::setMainSel(Step::selection* sel)
{
    curSel = sel;
    selection.update(sel);
}

void GuiWrapper::updateMainStep(bool draw_bonds, bool draw_cell)
{
    mainStep.update(curStep, draw_bonds, draw_cell);
}

void GuiWrapper::updateMainSelection()
{
    selection.update(curSel);
}

void GuiWrapper::addExtraData(GUI::Data* dat)
{
    extraData.insert(dat);
}

void GuiWrapper::delExtraData(GUI::Data* dat)
{
    extraData.erase(dat);
}

void GuiWrapper::updateViewUBO(void)
{
    if(rMatChanged){
        glBindBuffer(GL_UNIFORM_BUFFER, view_ubo);
        if(settings.perspective.val){
            glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(GUI::Mat),
                            (pMat*vMat*rMat).data());
        }else{
            glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(GUI::Mat),
                            (oMat*vMat*rMat).data());
        }
        glBufferSubData(GL_UNIFORM_BUFFER, sizeof(GUI::Mat), sizeof(GUI::Mat), rMat.data());
        rMatChanged = vMatChanged = pMatChanged = false;
    }else if (pMatChanged || vMatChanged){
        glBindBuffer(GL_UNIFORM_BUFFER, view_ubo);
        if(settings.perspective.val){
            glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(GUI::Mat),
                            (pMat*vMat*rMat).data());
        }else{
            glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(GUI::Mat),
                            (oMat*vMat*rMat).data());
        }
        vMatChanged = pMatChanged = false;
    }
}

void GuiWrapper::resizeViewMat(int w, int h)
{
    h==0?h=1:0;
    glViewport(0,0,w,h);
    float aspect = float(w)/h;
    oMat = guiMatMkOrtho(-10*aspect, 10*aspect, -10, 10, -100, 1000);
    pMat = guiMatMkPerspective(60, aspect, 0.001f, 1000);
    pMatChanged = true;
}

void GuiWrapper::zoomViewMat(int i)
{
    guiMatScale(vMat, i>0?1.1f:0.9f);
    vMatChanged = true;
}

void GuiWrapper::rotateViewMat(float x, float y, float z)
{
    guiMatRot(rMat, x, 0, 1, 0);
    guiMatRot(rMat, y, 1, 0, 0);
    guiMatRot(rMat, z, 0, 0, 1);
    rMatChanged = true;
}

void GuiWrapper::translateViewMat(float x, float y, float z)
{
    guiMatTranslate(vMat, x/10.f, y/10.f, z/10.f);
    vMatChanged = true;
}

void GuiWrapper::alignViewMat(alignDir d)
{
    switch (d) {
    case alignDir::x:
        rMat = {{ 0, 1, 0, 0,
                  0, 0, 1, 0,
                  1, 0, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::mx:
        rMat = {{ 0,-1, 0, 0,
                  0, 0, 1, 0,
                 -1, 0, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::y:
        rMat = {{-1, 0, 0, 0,
                  0, 0, 1, 0,
                  0, 1, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::my:
        rMat = {{ 1, 0, 0, 0,
                  0, 0, 1, 0,
                  0,-1, 0, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::z:
        rMat = {{ 1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1}};
        break;
    case alignDir::mz:
        rMat = {{-1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0,-1, 0,
                  0, 0, 0, 1}};
        break;
    }
    rMatChanged = true;
}

Mat GuiWrapper::getAxes()
{
    Mat tmp;
    tmp[0] =  Vec{rMat[0], rMat[1], rMat[2]};
    tmp[1] = -Vec{rMat[4], rMat[5], rMat[6]};
    tmp[2] =  Vec{rMat[8], rMat[9], rMat[10]};
    return tmp;
}
