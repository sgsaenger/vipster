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
        vbox = QVBoxLayout()
        vbox.addWidget(tooltip)
        vbox.addWidget(self.pickArea)
        self.setLayout(vbox)

    def setMol(self,mol):
        sel = mol.get_selection()
        if sel:
            br=0.52917721092
            at=[]
            vec = mol.get_vec()*mol.get_celldm()
            output =u'Atoms: '
            for j,i in enumerate(sel):
                at.append(mol.get_atom(i[0]))
                output+=str(i[0])+'('+at[j][0]+')'
                if j<len(sel)-1:
                    output+=', '
                else:
                    output+='\n'
            ids = [a[0] for a in sel]
            if len(sel)>1:
                #diff01 = sel[0][3]-sel[1][3]
                diff01 = at[0][1]+dot(vec,sel[0][1])-at[1][1]-dot(vec,sel[1][1])
                output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff01)*br,*ids[:2])
            if len(sel)>2:
                #3diff12 = sel[2][3]-sel[1][3]
                diff12 = at[1][1]+dot(vec,sel[1][1])-at[2][1]-dot(vec,sel[2][1])
                output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff12)*br,*ids[1:3])
                if len(sel)>3:
                    #diff23 = sel[2][3]-sel[3][3]
                    diff23 = at[2][1]+dot(vec,sel[2][1])-at[3][1]-dot(vec,sel[3][1])
                    output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff23)*br,*ids[2:])
                a012 = degrees(arccos(dot(diff01,-diff12)/(norm(diff01)*norm(diff12))))
                output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a012,*ids[:3])
            if len(sel)>3:
                a123 = degrees(arccos(dot(diff12,-diff23)/(norm(diff12)*norm(diff23))))
                output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a123,*ids[1:])
                c012 = cross(diff01,diff12)
                c123 = cross(diff12,diff23)
                if norm(c012)==0 or norm(c123)==0:
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0}\n'.format('not defined',*ids)
                else:
                    d0123 = degrees(arccos(dot(c012,c123)/(norm(c012)*norm(c123))))
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0:3.3f}°\n'.format(d0123,*ids)
            self.pickArea.setPlainText(output)
        else:
            self.pickArea.setPlainText('')
