#include <QMessageBox>
#include "atomlist.h"
#include "ui_atomlist.h"
#include "doubledelegate.h"
#include "vipsterapplication.h"

using namespace Vipster;

constexpr const char* inactiveFmt[] = { "Crystal", "Alat", "Angstrom", "Bohr" };
constexpr const char* activeFmt[] = {"Crystal (active)",
                                     "Alat (active)",
                                     "Angstrom (active)",
                                     "Bohr (active)"};

AtomList::AtomList(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::AtomList)
{
    ui->setupUi(this);

    // configure table
    ui->atomTable->setModel(&atomModel);
    connect(ui->atomTable->selectionModel(), &QItemSelectionModel::selectionChanged,
            this, &AtomList::selectionChanged);
    headerActions.push_back(new QAction{"Show type", ui->atomTable});
    headerActions.push_back(new QAction{"Show coordinates", ui->atomTable});
    headerActions.push_back(new QAction{"Show charges", ui->atomTable});
    headerActions.push_back(new QAction{"Show forces", ui->atomTable});
    headerActions.push_back(new QAction{"Show visibility", ui->atomTable});
    headerActions.push_back(new QAction{"Show constraints", ui->atomTable});
    for(auto& action: headerActions){
        action->setCheckable(true);
    }
    headerActions[0]->setChecked(true);
    headerActions[1]->setChecked(true);
    auto changeColumns = [&](){
        // triggered through changed columns
        int i=0;
        for(int j=0; j<headerActions.size(); ++j){
            i += headerActions[j]->isChecked() << j;
        }
        atomModel.setColumns(i);
    };
    for(auto& action: headerActions){
        connect(action, &QAction::toggled, this, changeColumns);
    }
    ui->atomTable->horizontalHeader()->setContextMenuPolicy(Qt::ActionsContextMenu);
    ui->atomTable->horizontalHeader()->addActions(headerActions);
    ui->atomTable->setItemDelegate(new DoubleDelegate{});

    // Connect ui elements
    connect(ui->displayButton, &QPushButton::toggled, ui->frame, &QWidget::setVisible);
    connect(ui->atomFmtButton, &QPushButton::clicked, this, &AtomList::fmtButtonHandler);
    connect(ui->atomFmtBox, &QComboBox::currentIndexChanged, this, &AtomList::fmtSelectionHandler);
    connect(ui->atomHelpButton, &QPushButton::clicked, [&](){
        QMessageBox::information(this, QString("About atoms"), Vipster::AtomsAbout);
    });

    // Connect to app state changes
    connect(&vApp, &Application::activeStepChanged, this, &AtomList::setActiveStep);
    connect(&vApp, &Application::stepChanged, this, &AtomList::updateStep);
    connect(&vApp, &Application::selChanged, this, &AtomList::updateSelection);
}

AtomList::~AtomList()
{
    delete ui;
}

void AtomList::setActiveStep(Step &step, Step::selection &sel)
{
    QSignalBlocker blockAtomFmt(ui->atomFmtBox);
    auto &fmtBox = *ui->atomFmtBox;

    // reset old fmt-string
    const auto oldFmt = fmtBox.currentIndex();
    fmtBox.setItemText(oldFmt, inactiveFmt[oldFmt]);

    // obtain formatter
    const auto fmt = step.getFmt();
    atomModel.setStep(step.asFmt(fmt));

    // mark fmt as active
    const auto ifmt = static_cast<int>(fmt)+2;
    fmtBox.setCurrentIndex(ifmt);
    fmtBox.setItemText(ifmt, activeFmt[ifmt]);

    // trigger further updates
    updateStep(step);
    updateSelection(sel);
}

void AtomList::updateStep(Step &step)
{
    // only update active step
    if (&step != &vApp.curStep()) return;

    atomModel.update();
}

void AtomList::updateSelection(Step::selection &sel)
{
    // only update selection of active step
    if (&sel != &vApp.curSel()) return;

    auto& table = *ui->atomTable;
    QSignalBlocker blockTable{table.selectionModel()};

    table.clearSelection();
    table.setSelectionMode(QAbstractItemView::MultiSelection);
    for(const auto& i: sel.getAtoms().indices){
        table.selectRow(static_cast<int>(i.first));
    }
    update(); // TODO: necessary?
    table.setSelectionMode(QAbstractItemView::ExtendedSelection);
}

void AtomList::selectionChanged()
{
    auto idx = ui->atomTable->selectionModel()->selectedRows();
    SelectionFilter filter{};
    filter.mode = SelectionFilter::Mode::Index;
    for(const auto& i: idx){
        filter.indices.emplace_back(static_cast<size_t>(i.row()), SizeVec{});
    }
    vApp.invokeOnSel([](Step::selection &sel, const SelectionFilter &filter){
        // TODO: sort out const-correctness
        sel = const_cast<Step&>(vApp.curStep()).select(filter);
    }, filter);
}

void AtomList::fmtSelectionHandler(int index)
{
    // TODO: sort out const-correctness
    atomModel.setStep(const_cast<Step&>(vApp.curStep()).asFmt(static_cast<AtomFmt>(index-2)));
}

void AtomList::fmtButtonHandler()
{
    QSignalBlocker blockFmtBox{ui->atomFmtBox};
    // reset old format string
    auto oldFmt = static_cast<int>(vApp.curStep().getFmt())+2;
    ui->atomFmtBox->setItemText(oldFmt, inactiveFmt[oldFmt]);

    // set new active format string
    auto ifmt = ui->atomFmtBox->currentIndex();
    ui->atomFmtBox->setItemText(ifmt, activeFmt[ifmt]);

    // modify actual Step
    auto fmt = static_cast<AtomFmt>(ifmt-2);
    vApp.invokeOnStep(&Step::setFmt, fmt, true);
    // TODO: sort out const-correctness
    ownStep = const_cast<Step&>(vApp.curStep()).asFmt(fmt);

    // reset selection
    SelectionFilter filter{};
    filter.mode = SelectionFilter::Mode::Index;
    filter.indices = vApp.curSel().getAtoms().indices;
    vApp.invokeOnSel([](Step::selection &sel, const SelectionFilter &filter){
        // TODO: sort out const-correctness
        sel = const_cast<Step&>(vApp.curStep()).select(filter);
    }, filter);
}
