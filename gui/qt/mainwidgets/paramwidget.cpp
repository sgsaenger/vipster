#include "paramwidget.h"
#include "ui_paramwidget.h"

#include <QMessageBox>
#include <QTreeWidget>
#include <QPlainTextEdit>
#include <QLineEdit>
#include <QAction>

using namespace Vipster;
using NameList = std::map<std::string, std::string>;

ParamWidget::ParamWidget(QWidget *parent) :
    BaseWidget(parent),
    ui(new Ui::ParamWidget)
{
    ui->setupUi(this);
    ui->valueArea->hide();
}

ParamWidget::~ParamWidget()
{
    delete ui;
}

void ParamWidget::clearParams()
{
    params.clear();
    ui->paramSel->clear();
}

void ParamWidget::registerParam(const std::string &name,
                                const Parameter &data)
{
    params.emplace_back(name, data);
    ui->paramSel->addItem(QString::fromStdString(
                          "(" +  data.getFmt()->command +
                           ") " + name
                         ));
    ui->paramSel->setCurrentIndex(ui->paramSel->count()-1);
}

void ParamWidget::on_paramSel_currentIndexChanged(int index)
{
    auto* layout = static_cast<QVBoxLayout*>(ui->valueArea->layout());
    // destroy old value-widgets
    while(auto *child = layout->takeAt(0)){
        if(child->layout()){
            while(auto *child_2 = child->layout()->takeAt(0)){
                layout->removeItem(child_2);
                child_2->widget()->deleteLater();
            }
            child->layout()->deleteLater();
        }else{
            child->widget()->deleteLater();
        }
    }
    // disable frame, show that nothing is loaded
    if(index<0){
        ui->noParam->show();
        ui->valueArea->hide();
        curParam = nullptr;
        return;
    }
    if(static_cast<size_t>(index) >= params.size()){
        throw Error("Invalid parameters selected");
    }
    // load values of new params
    ui->noParam->hide();
    ui->valueArea->show();
    curParam = &params.at(static_cast<size_t>(index)).second;
    // collect maps and vectors to combine complex widgets
    QVector<QString> vectors;
    QTreeWidget *tree{nullptr};
    for(auto& [name, pair]: *curParam){
        switch(pair.first.index()){
        case Parameter::i_strmap:
        {
            if(!tree){
                tree = setupTree();
                layout->insertWidget(0, tree);
            }
            QSignalBlocker block{tree};
            auto *top = new QTreeWidgetItem{{name.c_str()}};
            tree->addTopLevelItem(top);
            top->setToolTip(0, pair.second.c_str());
            for(const auto& [key, value]: std::get<NameList>(pair.first)){
                auto *child = new QTreeWidgetItem{{key.c_str(), value.c_str()}};
                child->setFlags(child->flags() | Qt::ItemIsEditable);
                top->addChild(child);
            }
            break;
        }
        case Parameter::i_strvec:
            vectors.push_back(QString::fromStdString(name));
            break;
        case Parameter::i_str:
        {
            // label + lineedit
            auto *row = new QHBoxLayout{};
            layout->addLayout(row);
            auto label = new QLabel{QString::fromStdString(name)+':'};
            row->addWidget(label);
            label->setToolTip(pair.second.c_str());
            auto edit = new QLineEdit{QString::fromStdString(
                        std::get<std::string>(pair.first))};
            row->addWidget(edit);
            edit->setToolTip(pair.second.c_str());
            connect(edit, &QLineEdit::editingFinished, [edit, &val=pair.first](){
                std::get<std::string>(val) = edit->text().toStdString();
            });
            break;
        }
        default:
            layout->addWidget(new QLabel{QString::fromStdString(name)+
                              " (unsupported setting)"});
            break;
        }
    }
    setupText(vectors);
}

QTreeWidget* ParamWidget::setupTree()
{
    auto *tree = new QTreeWidget{};
    tree->setContextMenuPolicy(Qt::ActionsContextMenu);
    tree->setColumnCount(2);
    tree->setHeaderLabels({"Parameter", "Value"});
    auto *addAction = new QAction{"Add entry", tree};
    tree->addAction(addAction);
    auto *delAction = new QAction{"Delete entry", tree};
    tree->addAction(delAction);
    // select entries
    connect(tree, &QTreeWidget::currentItemChanged,
            [&, delAction](QTreeWidgetItem *current, QTreeWidgetItem*){
        curItem = current;
        if(!curItem->parent()){
            delAction->setDisabled(true);
            curNL = &std::get<NameList>(curParam->at(curItem
                         ->text(0).toStdString()).first);
        }else{
            delAction->setEnabled(true);
            curNL = &std::get<NameList>(curParam->at(curItem->parent()
                         ->text(0).toStdString()).first);
            curKey = current->text(0).toStdString();
        }
    });
    // change entries
    connect(tree, &QTreeWidget::itemChanged, [&, tree](QTreeWidgetItem *item, int column){
        if(column == 0){
            // try to rename key
            auto newKey = item->text(column).toStdString();
            if(curNL->find(newKey) != curNL->end()){
                QMessageBox::warning(this, "Key already present",
                    QString{"Cannot rename key \""}+curKey.c_str()+"\" to \""+
                    newKey.c_str()+"\" because the latter is already present.");
                QSignalBlocker block{tree};
                item->setText(0, QString::fromStdString(curKey));
            }else{
                auto node = curNL->extract(curKey);
                node.key() = newKey;
                curNL->insert(std::move(node));
            }
        }else{
            // change value
            curNL->at(item->text(0).toStdString()) = item->text(1).toStdString();
        }
    });
    // add new entries
    connect(addAction, &QAction::triggered, [&](){
        std::string key = "newKey";
        if(curNL->find(key) != curNL->end()){
            uint8_t idx{0};
            std::string key2 = key + "0";
            while(curNL->find(key2) != curNL->end()){
                key2 = key + std::to_string(++idx);
            }
            key = key2;
        }
        auto *child = new QTreeWidgetItem{{key.c_str(), (*curNL)[key].c_str()}};
        child->setFlags(child->flags() | Qt::ItemIsEditable);
        if(!curItem->parent()){
            curItem->addChild(child);
        }else{
            curItem->parent()->addChild(child);
        }
    });
    // delete Entries
    connect(delAction, &QAction::triggered, [&](){
        curNL->erase(curNL->find(curItem->text(0).toStdString()));
        delete curItem;
    });
    return tree;
}

void ParamWidget::setupText(const QVector<QString> &vectors)
{
    // if vectors are present, add a textedit,
    // and a combobox if multiple vectors are needed
    auto* layout = static_cast<QVBoxLayout*>(ui->valueArea->layout());
    QPlainTextEdit *text = vectors.empty() ? nullptr : new QPlainTextEdit{};
    QComboBox *vecSel = vectors.size()>1 ? new QComboBox{} : nullptr;
    if(text){
        layout->insertWidget(0, text);
        connect(text, &QPlainTextEdit::textChanged, [&, text](){
        auto pos = curParam->find(curVec);
        if (pos == curParam->end()) return;
        auto &value = std::get<std::vector<std::string>>(pos->second.first);
        value.clear();
        auto tmp = text->toPlainText();
        if(tmp.size()){
            for(const auto &line: tmp.split('\n')){
                value.push_back(line.toStdString());
            }
        }
    });
    }
    if(vecSel){
        vecSel->addItems(QStringList::fromVector(vectors));
        layout->insertWidget(0, vecSel);
        auto fillText = [&, text, vecSel](int i){//const QString& s){
            QSignalBlocker block{text};
            curVec = vecSel->currentText().toStdString();
            auto &[value, doc] = curParam->at(curVec);
            text->clear();
            QStringList list{};
            const auto &vec = std::get<std::vector<std::string>>(value);
            for(const auto& line: vec){
                list.append(QString::fromStdString(line));
            }
            text->setPlainText(list.join('\n'));
            text->setToolTip(doc.c_str());
        };
        fillText(0);
        connect(vecSel, QOverload<int>::of(&QComboBox::currentIndexChanged), fillText);
    }else if(text){
        curVec = vectors[0].toStdString();
        auto label = new QLabel{QString::fromStdString(curVec)+':'};
        layout->insertWidget(0, label);
        label->setBuddy(text);
        label->setToolTip(curParam->at(curVec).second.c_str());
        QStringList list{};
        for(const auto& line: std::get<std::vector<std::string>>(curParam->at(curVec).first)){
            list.append(QString::fromStdString(line));
        }
        text->setPlainText(list.join('\n'));
        text->setToolTip(curParam->at(curVec).second.c_str());
    }
}

void ParamWidget::on_pushButton_clicked()
{
    QMessageBox::information(this, QString("About parameter presets"),
                             Vipster::ParametersAbout);
}
