# -*- coding: utf-8 -*-
from ..molecule import Molecule
from ..settings import config

from re import split as regex
from numpy import cos,isclose,diag

name = 'CPMD Input File'
extension = 'cpi'
argument = 'cpmd'

param = None
writer = None

def parser(name,data):
    """ Parse CPMD Input file """
    tmol = Molecule(name)
    i=0
    ibrav='14'
    symmetries={'ISOLATED':'0',
                'CUBIC':'1',
                'FACE CENTERED CUBIC':'2',
                'FCC':'2',
                'BODY CENTERED CUBIC':'3',
                'BCC':'3',
                'HEXAGONAL':'4',
                'TRIGONAL':'5',
                'RHOMBOHEDRAL':'5',
                'TETRAGONAL':'6',
                'BODY CENTERED TETRAGONAL':'7',
                'BCT':'7',
                'ORTHOROMBIC':'8',
                'MONOCLINIC':'12',
                'TRICLINIC':'14'}
    scale=None
    unitvec=[[1,0,0],[0,1,0],[0,0,1]]
    scalevec=[[1,0,0],[0,1,0],[0,0,1]]
    angstrom=False
    while i < len(data):
        if data[i][0]=='!' or data[i].split()==[]:
            i+=1
        elif "&SYSTEM" in data[i]:
            #parse cell geometry!
            while not "&END" in data[i]:
                if "ANGSTROM" in data[i]:
                    scale="angstrom"
                elif "SCALE" in data[i]:
                    if "CARTESIAN" in data[i]:
                        scale="alat"
                    else:
                        scale="crystal"
                    for token in data[i].strip().split():
                        if 'S=' in token:
                            scalevec[0][0]=1/float(token.strip('S='))
                            scalevec[1][1]=1/float(token.strip('S='))
                            scalevec[2][2]=1/float(token.strip('S='))
                        elif 'SX=' in token:
                            scalevec[0][0]=1/float(token.strip('SX='))
                        elif 'SY=' in token:
                            scalevec[1][1]=1/float(token.strip('SY='))
                        elif 'SZ=' in token:
                            scalevec[2][2]=1/float(token.strip('SZ='))
                elif "KPOINTS" in data[i]:
                    if "MONKHORST" in data[i]:
                        kpoints=data[i+1].split()
                        if "SHIFT=" in data[i]:
                            kpoints+=data[i].split("SHIFT=")[1].split()[0:3]
                        else:
                            kpoints+=['0','0','0']
                        tmol.setKpoints('active','mpg')
                        tmol.setKpoints('mpg',kpoints)
                        i+=1
                    else:
                        kpoints=[]
                        if "SCALED" in data[i]:
                            crystal=True
                        if "BANDS" in data[i]:
                            bands=True
                            j=0
                            line=data[i+j].split()
                            while not all([float(i)==0 for i in line]):
                                kpoints.append(line)
                                j+=1
                                line=data[i+j].split()
                            i+=j
                        else:
                            nk = data[i+1]
                            for j in range(nk):
                                kpoints.append(data[i+j+1].split())
                            i+=nk
                elif "SYMMETRY" in data[i]:
                    arg = data[i+1].strip()
                    if arg.isdigit():
                        ibrav=arg
                    else:
                        ibrav=symmetries[arg]
                    i+=1
                elif "CELL" in data[i]:
                    cell=list(map(float,data[i+1].split()))
                    if "VECTORS" in data[i]:
                        ibrav='-2'
                    else:
                        #cell = [a,b/a,c/a,cos(alpha),cos(beta),cos(gamma)]
                        a,b,c,cosal,cosbeta,cosgam=cell
                        tmol.setCellDim(a)
                        if "ABSOLUTE" in data[i]:
                            #cell = [a,b,c,...]
                            b=b/a
                            c=c/a
                        elif "DEGREE" in data[i]:
                            #cell = [...,alpha,beta,gamma]
                            cosal,cosbeta,cosgam=map(cos,[cosal,cosbeta,cosgam])
                    i+=1
                i+=1
        elif '&ATOMS' in data[i]:
            #parse atom coordinates
            j=i+1
            while not "&END" in data[j]:
                if data[j][0]=='*':
                    atype=regex('[-_.]',data[j][1:].strip())[0]
                    nat=int(data[j+2])
                    tmol.newAtoms(nat)
                    for k in range(nat):
                        tmol.setAtom(k,atype,data[j+3+k].split())
                    j+=2+nat
                j+=1
                pass
            i+=j
            pass
        else:
            i+=1
    #apply cell vec (scale SX etc.)
    if scalevec!=unitvec:
        tmol.setVec(scalevec,scale=True)
        tmol.setVec(unitvec)
    #parse cell parameters
    if ibrav=='-2':
        tmol.setVec([cell[0:3],cell[3:6],cell[6:9]])
    if ibrav == '1':
        #simple cubic
        pass
    elif ibrav == '2':
        #face centered cubic
        tmol.setVec([[-0.5,0,0.5],[0,0.5,0.5],[-0.5,0.5,0]])
    elif ibrav == '3':
        #body centered cubic
        tmol.setVec([[0.5,0.5,0.5],[-0.5,0.5,0.5],[-0.5,-0.5,0.5]])
    elif ibrav == '4':
        #hexagonal
        tmol.setVec([[1,0,0],[-0.5,sqrt(3)*0.5,0],[0,0,c]])
    elif ibrav == '5':
        #trigonal
        tx=sqrt((1-alpha)/2)
        ty=sqrt((1-alpha)/6)
        tz=sqrt((1+2*alpha)/3)
        tmol.setVec([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]])
    elif ibrav == '6':
        #simple tetragonal
        tmol.setVec([[1,0,0],[0,1,0],[0,0,c]])
    elif ibrav == '7':
        #body centered tetragonal
        tmol.setVec([[0.5,-0.5,c*0.5],[0.5,0.5,c*0.5],[-0.5,-0.5,c*0.5]])
    elif ibrav=='0' or ibrav=='8':
        #simple orthorhombic
        tmol.setVec([[1,0,0],[0,b,0],[0,0,c]])
    elif ibrav == '12':
        #simple monoclinic
        tmol.setVec([[1,0,0],[b*alpha,b*sqrt(1-alpha**2),0],[0,0,c]])
    elif ibrav == '14':
        #triclinic
        singam = sqrt(1-gamma**2)
        tmol.setVec([[1,0,0],[b*gamma,b*singam,0],
                [c*beta,c*(alpha-beta*gamma)/singam,c*sqrt(1+2*alpha*beta*gamma-alpha*alpha-beta*beta-gamma*gamma)/singam]])
    #scale atoms/cell
    if scale:
        tmol.scaleAtoms(scale)
        tmol.setCellDim(tmol.getCellDim(),fmt=scale)
    return tmol,None
