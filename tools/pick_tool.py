#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from numpy.linalg import norm
from numpy import degrees,arccos,dot,cross

####################################
# Selected Atoms
####################################

class Picker(QWidget):

    def __init__(self,parent):
        super(Picker,self).__init__()

        tooltip = QLabel('Pick up to 4 atoms:')
        self.pickArea = QTextEdit()
        self.pickArea.setReadOnly(True)
        self.pickWarn = QLabel()
        vbox = QVBoxLayout()
        vbox.addWidget(tooltip)
        vbox.addWidget(self.pickArea)
        vbox.addWidget(self.pickWarn)
        self.setLayout(vbox)

    def setSel(self,sel):
        br=0.52917721092
        self.pickWarn.setText('')
        if len(sel)==0:
            self.pickArea.setPlainText('')
        else:
            output ='Atoms: '+str([a[1] for a in sel])+'\n'
            output+='Types: '+str([a[2] for a in sel])+'\n'
            ids = [a[1] for a in sel]
            if len(sel)>1:
                diff01 = sel[0][3]-sel[1][3]
                output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff01)*br,*ids[:2])
            if len(sel)>2:
                diff12 = sel[2][3]-sel[1][3]
                output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff12)*br,*ids[1:3])
                if len(sel)>3:
                    diff23 = sel[2][3]-sel[3][3]
                    output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff23)*br,*ids[2:])
                a012 = degrees(arccos(dot(diff01,diff12)/(norm(diff01)*norm(diff12))))
                output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a012,*ids[:3])
            if len(sel)>3:
                a123 = degrees(arccos(dot(diff12,diff23)/(norm(diff12)*norm(diff23))))
                output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a123,*ids[1:])
                c012 = cross(diff01,diff12)
                c123 = cross(diff12,diff23)
                if norm(c012)==0 or norm(c123)==0:
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0}\n'.format('not defined',*ids)
                else:
                    d0123 = degrees(arccos(dot(c012,c123)/(norm(c012)*norm(c123))))
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0:3.3f}°\n'.format(d0123,*ids)
            self.pickArea.setPlainText(output)

    def setMol(self,mol):
        if self.pickArea.toPlainText():
            self.pickWarn.setText('Data has changed!')
