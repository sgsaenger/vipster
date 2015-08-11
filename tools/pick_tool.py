#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from PyQt4.QtGui import *
from numpy.linalg import norm
from numpy import degrees,arccos,dot,cross

####################################
# Selected Atoms
####################################

class Picker(QTextEdit):

    def __init__(self,parent):
        super(Picker,self).__init__()

        self.setReadOnly(True)

    def setMol(self,mol):
        sel = mol.get_selection()
        if sel:
            at=[]
            vec = mol.get_vec()*mol.get_celldm('angstrom')
            output =u'Atoms: '
            for j,i in enumerate(sel):
                at.append(mol.get_atom(i[0],'angstrom'))
                output+=str(i[0])+'('+at[j][0]+')'
                if j<len(sel)-1:
                    output+=', '
                else:
                    output+='\n'
            ids = [a[0] for a in sel]
            if len(sel)>1:
                diff01 = at[0][1]+dot(sel[0][1],vec)-at[1][1]-dot(sel[1][1],vec)
                output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff01),*ids[:2])
            if len(sel)>2:
                diff12 = at[1][1]+dot(sel[1][1],vec)-at[2][1]-dot(sel[2][1],vec)
                output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff12),*ids[1:3])
                if len(sel)>3:
                    diff23 = at[2][1]+dot(sel[2][1],vec)-at[3][1]-dot(sel[3][1],vec)
                    output+=u'Dist {1}-{2}: {0:3.3f} Å\n'.format(norm(diff23),*ids[2:])
                if norm(diff01)>0 and norm(diff12)>0:
                    a012 = degrees(arccos(dot(diff01,-diff12)/(norm(diff01)*norm(diff12))))
                    output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a012,*ids[:3])
                else:
                    output+=u'Angle {1}-{2}-{3}: {0}°\n'.format('not defined',*ids[:3])
            if len(sel)>3:
                if norm(diff12)>0 and norm(diff23)>0:
                    a123 = degrees(arccos(dot(diff12,-diff23)/(norm(diff12)*norm(diff23))))
                    output+=u'Angle {1}-{2}-{3}: {0:3.3f}°\n'.format(a123,*ids[1:])
                else:
                    output+=u'Angle {1}-{2}-{3}: {0}°\n'.format('not defined',*ids[1:])
                c012 = cross(diff01,diff12)
                c123 = cross(diff12,diff23)
                if norm(c012)>0 or norm(c123)>0:
                    d0123 = degrees(arccos(dot(c012,c123)/(norm(c012)*norm(c123))))
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0:3.3f}°\n'.format(d0123,*ids)
                else:
                    output+=u'Dihedral {1}-{2}-{3}-{4}: {0}\n'.format('not defined',*ids)
            #self.pickArea.setPlainText(output)
            self.setPlainText(output)
        else:
            #self.pickArea.setPlainText('')
            self.setPlainText('')
