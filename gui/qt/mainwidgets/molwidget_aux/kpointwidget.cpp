#include "kpointwidget.h"
#include "ui_kpointwidget.h"
#include "vipsterapplication.h"

using namespace Vipster;

constexpr const char* fmtNames[] = {"Gamma", "Monkhorst-Pack grid", "Discrete"};

KPointWidget::KPointWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::KPointWidget)
{
    ui->setupUi(this);

    // Connect to app state changes
    connect(&vApp, &Application::activeMolChanged, this, &KPointWidget::setActiveMol);
    connect(&vApp, &Application::molChanged, this, &KPointWidget::updateMol);

    // Connect ui elements
    connect(ui->displayButton, &QPushButton::toggled, ui->frame, &QWidget::setVisible);
    ui->frame->setVisible(ui->displayButton->isChecked());
    connect(ui->kFmtButton, &QPushButton::clicked, this, [this](){
        vApp.invokeOnMol([](Molecule &m, KPoints::Fmt fmt){
            m.kpoints.active = fmt;
        }, static_cast<KPoints::Fmt>(ui->activeKpoint->currentIndex()));
    });
    connect(ui->activeKpoint, &QComboBox::currentIndexChanged, ui->kpointStack, &QStackedWidget::setCurrentIndex);

    // Monkhorst-Pack elements
    connect(ui->mpg_x, &QSpinBox::valueChanged,
            [](int i){vApp.invokeOnMol([](Molecule &m, int i){m.kpoints.mpg.x = i;}, i);});
    connect(ui->mpg_y, &QSpinBox::valueChanged,
            [](int i){vApp.invokeOnMol([](Molecule &m, int i){m.kpoints.mpg.y = i;}, i);});
    connect(ui->mpg_z, &QSpinBox::valueChanged,
            [](int i){vApp.invokeOnMol([](Molecule &m, int i){m.kpoints.mpg.z = i;}, i);});
    connect(ui->mpg_x_off, &QDoubleSpinBox::valueChanged,
            [](double d){vApp.invokeOnMol([](Molecule &m, double d){m.kpoints.mpg.sx = d;}, d);});
    connect(ui->mpg_y_off, &QDoubleSpinBox::valueChanged,
            [](double d){vApp.invokeOnMol([](Molecule &m, double d){m.kpoints.mpg.sy = d;}, d);});
    connect(ui->mpg_z_off, &QDoubleSpinBox::valueChanged,
            [](double d){vApp.invokeOnMol([](Molecule &m, double d){m.kpoints.mpg.sz = d;}, d);});

    // Discrete kPoint table
    // TODO: add tool tips to explain modes
    auto newAction = new QAction{"New K-Point"};
    ui->discretetable->addAction(newAction);
    connect(newAction, &QAction::triggered,
            [](){ vApp.invokeOnMol( [](Molecule &m){
                    m.kpoints.discrete.kpoints.emplace_back();
            });
    });
    auto delAction = new QAction{"Delete K-Point"};
    ui->discretetable->addAction(delAction);
    delAction->setDisabled(true);
    connect(delAction, &QAction::triggered,
            this,
            [this](){
                auto sel = ui->discretetable->selectedItems();
                if (!sel.empty()) {
                    vApp.invokeOnMol( [](Molecule &m, int i){
                        auto &kpoints = m.kpoints.discrete.kpoints;
                        kpoints.erase(kpoints.begin() + i);
                    }, sel[0]->row());
                }
            });
    connect(ui->discretetable, &QTableWidget::itemSelectionChanged, delAction,
            [delAction, this](){
                auto sel = ui->discretetable->selectedItems();
                delAction->setEnabled(!sel.empty());
            });
    connect(ui->discretetable, &QTableWidget::cellChanged,
            [this](int row, int col){
                auto *cell = ui->discretetable->item(row, col);
                vApp.invokeOnMol([](Molecule &m, int row, int col, double val){
                    if (col < 3) {
                        m.kpoints.discrete.kpoints[row].pos[col] = val;
                    } else {
                        m.kpoints.discrete.kpoints[row].weight = val;
                    }
                }, row, col, cell->text().toDouble());
            });
    // TODO: support multiple modes
    connect(ui->bands, &QCheckBox::stateChanged,
            [](int i){vApp.invokeOnMol([](Molecule &m, bool b){
                if (b) {
                    m.kpoints.discrete.properties |= KPoints::Discrete::band;
                } else {
                    m.kpoints.discrete.properties ^= KPoints::Discrete::band;
                }
            }, i);});
    connect(ui->crystal, &QCheckBox::stateChanged,
            [](int i){vApp.invokeOnMol([](Molecule &m, bool b){
                if (b) {
                    m.kpoints.discrete.properties |= KPoints::Discrete::crystal;
                } else {
                    m.kpoints.discrete.properties ^= KPoints::Discrete::crystal;
                }
            }, i);});
}

KPointWidget::~KPointWidget()
{
    delete ui;
}

void KPointWidget::setActiveMol(const Molecule &m)
{
    // show active format
    ui->activeKpoint->setCurrentIndex(static_cast<int>(m.kpoints.active));

    // Update other content
    updateMol(m);
}

void KPointWidget::updateMol(const Molecule &m)
{
    // Only act on active molecule
    if (&m != &vApp.curMol()) return;

    // signal active format
    for (int i=0; i<3; ++i) {
        if (i == static_cast<int>(m.kpoints.active)) {
            ui->activeKpoint->setItemText(i, QString{fmtNames[i]} + " (active)");
        } else {
            ui->activeKpoint->setItemText(i, fmtNames[i]);
        }
    }

    // fill MPG infos
    ui->mpg_x->setValue(m.kpoints.mpg.x);
    ui->mpg_y->setValue(m.kpoints.mpg.y);
    ui->mpg_z->setValue(m.kpoints.mpg.z);
    ui->mpg_x_off->setValue(m.kpoints.mpg.sx);
    ui->mpg_y_off->setValue(m.kpoints.mpg.sy);
    ui->mpg_z_off->setValue(m.kpoints.mpg.sz);

    // fill discrete table
    ui->crystal->setChecked(m.kpoints.discrete.properties & KPoints::Discrete::crystal);
    ui->bands->setChecked(m.kpoints.discrete.properties & KPoints::Discrete::band);
    QSignalBlocker blockTable{ui->discretetable};
    auto &table = *ui->discretetable;
    const auto &points = m.kpoints.discrete.kpoints;
    auto count = static_cast<int>(points.size());
    table.clear();
    table.setRowCount(count);
    for (int i=0; i<count; ++i) {
        const auto &kp = points[i];
        for (int j=0; j<3; ++j) {
            table.setItem(i, j, new QTableWidgetItem(QString::number(kp.pos[j])));
        }
        table.setItem(i, 3, new QTableWidgetItem(QString::number(kp.weight)));
    }
}
