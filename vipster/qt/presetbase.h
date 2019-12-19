#ifndef CONFIGBASE_H
#define CONFIGBASE_H

#include <QWidget>
#include "fileio.h"

class PresetBase: public QWidget
{
    Q_OBJECT

public:
    explicit PresetBase(QWidget *parent = nullptr);
    virtual ~PresetBase() = default;
    virtual void setPreset(Vipster::IO::Preset *c)=0;
};

#endif // CONFIGBASE_H
