#include "pickwidget.h"
#include "ui_pickwidget.h"
using namespace Vipster;

PickWidget::PickWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::PickWidget)
{
    ui->setupUi(this);
}

PickWidget::~PickWidget()
{
    delete ui;
}

inline void printDist(QPlainTextEdit& text,
                      QString idx1, QString idx2, float dist)
{
    QString tmp = "Dist " + idx1 + '-' + idx2 + ": " +
            QString::number(static_cast<double>(dist)) + " Å";
    text.appendPlainText(tmp);
}

inline void printAngle(QPlainTextEdit& text,
                       QString idx1, QString idx2, QString idx3, float ang)
{
    QString tmp = "Angle " + idx1 + '-' +
            idx2 + '-' + idx3 + ": " +
            QString::number(static_cast<double>(ang)) + "°";
    text.appendPlainText(tmp);
}

inline void printDihed(QPlainTextEdit& text,
                       QString idx1, QString idx2, QString idx3, QString idx4, float dihed)
{
    QString tmp = "Dihed " +
            idx1 + '-' + idx2 + '-' +
            idx3 + '-' + idx4 + ": " +
            QString::number(static_cast<double>(dihed)) + "°";
    text.appendPlainText(tmp);
}

void PickWidget::updateWidget(guiChange_t change)
{
    if((change & (GuiChange::atoms|GuiChange::cell|GuiChange::selection)) == 0u){
        return;
    }
    const auto& curSel = master->curSel->asFmt(AtomFmt::Angstrom);
    auto& text = *ui->PickText;
    text.setPlainText("Atoms:");
    auto it = curSel.begin();
    size_t nat = 0;
    std::vector<QString> names;
    std::vector<Vec> coords;
    auto fmt = curSel.getFormatter(AtomFmt::Crystal, AtomFmt::Angstrom);
    while(it != curSel.end()){
        int count = 0;
        const auto& pair = it.getFilterPair();
        for(const auto& off: pair.second){
            names.push_back(QString::number(pair.first)+QString{count,'\''});
            if(nat<4){
                // save four coordinates in unwrapped form
                coords.push_back(it->coord + fmt(Vec{(float)off[0],(float)off[1],(float)off[2]}));
            }
            if(off != SizeVec{0,0,0}){
                text.appendPlainText(names.back()+'('+
                                     QString::fromStdString(it->name)+") <"+
                                     QString::number(off[0])+','+
                                     QString::number(off[1])+','+
                                     QString::number(off[2])+'>');
            }else{
                text.appendPlainText(names.back()+'('+
                                     QString::fromStdString(it->name)+')');
            }
            count++;
            nat++;
        }
        ++it;
    }
    if(nat>4){
        //Don't display additional information when too many atoms are selected
        return;
    }
    /*
     *  |-----|------|-----|-------|
     *  | Nat | Dist | Ang | Dihed |
     *  |-----|------|-----|-------|
     *  | 2   | 1    | -   | -     |
     *  |-----|------|-----|-------|
     *  | 3   | 3    | 3   | -     |
     *  |-----|------|-----|-------|
     *  | 4   | 7    | -   | 6     |
     *  |-----|------|-----|-------|
     */
    if(nat>1){
        auto diff01 = coords[0] - coords[1];
        auto dist01 = Vec_length(diff01);
        printDist(text, names[0], names[1], dist01);
        if(nat>2){
            auto diff02 = coords[0] - coords[2];
            auto diff12 = coords[1] - coords[2];
            auto dist02 = Vec_length(diff02);
            auto dist12 = Vec_length(diff12);
            printDist(text, names[0], names[2], dist02);
            printDist(text, names[1], names[2], dist12);
            if(nat == 3){
                auto ang012 = (std::acos(Vec_dot(diff01, diff12) / (dist01 * dist12)))
                              * rad2deg;
                auto ang120 = (std::acos(Vec_dot(diff02, diff12) / (dist02 * dist12)))
                              * rad2deg;
                auto ang201 = (std::acos(Vec_dot(diff01, diff02) / (dist01 * dist02)))
                              * rad2deg;
                printAngle(text, names[0], names[1], names[2], ang012);
                printAngle(text, names[1], names[2], names[0], ang120);
                printAngle(text, names[2], names[0], names[1], ang201);
            }else{
                auto diff03 = coords[0] - coords[3];
                auto diff13 = coords[1] - coords[3];
                auto diff23 = coords[2] - coords[3];
                auto dist03 = Vec_length(diff03);
                auto dist13 = Vec_length(diff13);
                auto dist23 = Vec_length(diff23);
                printDist(text, names[0], names[3], dist03);
                printDist(text, names[1], names[3], dist13);
                printDist(text, names[2], names[3], dist23);
                auto cross012 = Vec_cross(diff01, diff12);
                auto cross013 = Vec_cross(diff01, diff13);
                auto cross023 = Vec_cross(diff02, diff23);
                auto cross123 = Vec_cross(diff12, diff23);
                auto len012 = Vec_length(cross012);
                auto len013 = Vec_length(cross013);
                auto len023 = Vec_length(cross023);
                auto len123 = Vec_length(cross123);
                auto dihed0123 = (std::acos(Vec_dot(cross012, cross123) / (len012 * len123)))
                                 * rad2deg;
                auto dihed0132 = (std::acos(Vec_dot(cross013, cross123) / (len013 * len123)))
                                 * rad2deg;
                auto dihed1023 = (std::acos(Vec_dot(cross012, cross023) / (len012 * len023)))
                                 * rad2deg;
                auto dihed1230 = (std::acos(Vec_dot(cross023, cross123) / (len023 * len123)))
                                 * rad2deg;
                auto dihed1032 = (std::acos(Vec_dot(cross013, cross023) / (len013 * len023)))
                                 * rad2deg;
                auto dihed2013 = (std::acos(Vec_dot(cross012, cross013) / (len012 * len013)))
                                 * rad2deg;
                printDihed(text, names[0], names[1], names[2], names[3], dihed0123);
                printDihed(text, names[0], names[1], names[3], names[2], dihed0132);
                printDihed(text, names[1], names[0], names[2], names[3], dihed1023);
                printDihed(text, names[1], names[2], names[3], names[0], dihed1230);
                printDihed(text, names[1], names[0], names[3], names[2], dihed1032);
                printDihed(text, names[2], names[0], names[1], names[3], dihed2013);
            }
        }
    }
}
