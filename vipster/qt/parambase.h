#ifndef PARAMBASE_H
#define PARAMBASE_H

#include <QWidget>
#include "fileio.h"

class ParamBase: public QWidget
{
    Q_OBJECT

public:
    explicit ParamBase(QWidget *parent = nullptr);
    virtual ~ParamBase() = default;
    virtual void setParam(Vipster::IO::Parameter *p)=0;
};

#endif // PARAMBASE_H
