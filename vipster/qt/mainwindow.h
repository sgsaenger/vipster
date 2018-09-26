#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QAbstractButton>
#include <QDir>
#include <QTimer>
#include <vector>
#include "io.h"
#include "stepsel.h"
#include "../common/guiwrapper.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QString path, QWidget *parent = nullptr);
    explicit MainWindow(QString path,
                        std::vector<Vipster::IO::Data> &&d,
                        QWidget *parent = nullptr);
    ~MainWindow() override;
    Vipster::Molecule* curMol{nullptr};
    Vipster::StepProper* curStep{nullptr};
    Vipster::StepSelection* curSel{nullptr};
    Vipster::AtomFmt getFmt();
    void setFmt(int i, bool apply, bool scale);
    void updateWidgets(uint8_t change);
    void newData(Vipster::IO::Data&& d);
    std::list<Vipster::Molecule> molecules;
    std::list<std::unique_ptr<const Vipster::BaseData>> data;
    const std::vector<std::pair<Vipster::IOFmt, std::unique_ptr<Vipster::IO::BaseParam>>>& getParams() const noexcept;
    const std::vector<std::pair<Vipster::IOFmt, std::unique_ptr<Vipster::IO::BaseConfig>>>& getConfigs() const noexcept;
    void addExtraData(Vipster::GUI::Data* dat);
    void delExtraData(Vipster::GUI::Data* dat);
    const Vipster::GUI::GlobalData& getGLGlobals();

public slots:
    void setMol(int i);
    void setStep(int i);
    void stepBut(QAbstractButton *but);
    void about(void);
    void newMol();
    void loadMol();
    void saveMol();
    void editAtoms(void);
    void loadParam();
    void saveParam();
    void loadConfig();
    void saveConfig();
    void saveScreenshot();

private:
    void setupUI(void);

    Vipster::AtomFmt fmt{Vipster::AtomFmt::Bohr};
    Ui::MainWindow *ui;
    QDir path{};
    QTimer playTimer{};
    std::vector<QDockWidget*> baseWidgets;
    std::vector<QDockWidget*> toolWidgets;
};

class BaseWidget: public QWidget{
    Q_OBJECT

public:
    BaseWidget(QWidget *parent = nullptr);
    void triggerUpdate(uint8_t change);
    virtual void updateWidget(uint8_t){}
    virtual ~BaseWidget()=default;

protected:
    bool updateTriggered{false};
    MainWindow* master;
};

#endif // MAINWINDOW_H
