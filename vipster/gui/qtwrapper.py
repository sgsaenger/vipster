# -*- coding: utf-8 -*-

try:
    from PyQt5.QtGui import *
    from PyQt5.QtWidgets import *
    from PyQt5.QtOpenGL import *
    from PyQt5.QtCore import Qt,QRectF,QTimer,QSettings,QByteArray,qVersion
    QGLShader = QOpenGLShader
    QGLShaderProgram = QOpenGLShaderProgram
    try:
        from PyQt5.QtCore import QString
    except:
        QString = str
except:
    from PyQt4.QtGui import *
    from PyQt4.QtOpenGL import *
    from PyQt4.QtCore import Qt,QRectF,QTimer,QSettings,QByteArray,qVersion
    try:
        from PyQt4.QtCore import QString
    except:
        QString = str
