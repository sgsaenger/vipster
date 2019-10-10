#ifndef CONFIGBASE_H
#define CONFIGBASE_H

#include <QWidget>
#include "io.h"

class ConfigBase: public QWidget
{
    Q_OBJECT

public:
    explicit ConfigBase(QWidget *parent = nullptr);
    virtual ~ConfigBase() = default;
    virtual void setConfig(Vipster::IO::BaseConfig *c)=0;
};

#endif // CONFIGBASE_H
