# -*- coding: utf-8 -*-
from ..molecule import Molecule
from ..settings import config

from re import split as regex
from numpy import cos,isclose,diag
from collections import OrderedDict
from copy import deepcopy

name = 'CPMD Input'
extension = 'cpi'
argument = 'cpmd'

param = {"default":OrderedDict([("type","cpmd"),
                    ("name","New CPMD parameters"),
                    ("&INFO",""),
                    ("&CPMD",""),
                    ("&SYSTEM",""),
                    ("&PIMD",""),
                    ("&PATH",""),
                    ("&ATOMS",""),
                    ("&DFT",""),
                    ("&PROP",""),
                    ("&BASIS",""),
                    ("&RESP",""),
                    ("&PTDDFT",""),
                    ("&LINRES",""),
                    ("&HARDNESS",""),
                    ("&TDDFT",""),
                    ("&QMMMM",""),
                    ("&CLAS",""),
                    ("&EXTE",""),
                    ("&VDW","")])}

def parser(name,data):
    """ Parse CPMD Input file """
    tmol = Molecule(name)
    tparam = deepcopy(param["default"])
    tparam["name"]=name
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
        elif '&' in data[i] and data[i].strip()!="&END":
            #start recording namelist
            nl=data[i].strip()
            i+=1
            while not "&END" in data[i]:
                if nl=="&SYSTEM" and "ANGSTROM" in data[i]:
                    scale="angstrom"
                elif nl=="&SYSTEM" and "SCALE" in data[i]:
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
                elif nl=="&SYSTEM" and "KPOINTS" in data[i]:
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
                        crystal=False
                        bands=False
                        if "SCALED" in data[i]:
                            crystal=True
                        if "BANDS" in data[i]:
                            bands=True
                            j=1
                            line=data[i+j].split()
                            while not all([float(k)==0 for k in line]):
                                kpoints.append(line[1:4]+[line[0]])
                                kpoints.append(line[4:]+['0'])
                                j+=1
                                line=data[i+j].split()
                            i+=j
                        else:
                            nk = int(data[i+1])
                            for j in range(nk):
                                kpoints.append(data[i+j+2].split())
                            i+=nk
                        tmol.setKpoints('active','discrete')
                        tmol.setKpoints('discrete',kpoints)
                        tmol.setKpoints('options',{'crystal':crystal,'bands':bands})
                elif nl=="&SYSTEM" and "SYMMETRY" in data[i]:
                    arg = data[i+1].strip()
                    if arg.isdigit():
                        ibrav=arg
                    else:
                        ibrav=symmetries[arg]
                    i+=1
                elif nl=="&SYSTEM" and "CELL" in data[i]:
                    cell=list(map(float,data[i+1].split()))
                    if "VECTORS" in data[i]:
                        while len(cell)<9:
                            i+=1
                            cell+=list(map(float,data[i+1].split()))
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
                elif nl=="&ATOMS" and data[i][0]=='*':
                    atype=regex('[-_.]',data[i][1:].strip())[0]
                    tmol.pse[atype]['CPPP']=data[i][1:].strip()
                    tmol.pse[atype]['CPNL']=data[i+1].strip()
                    nat=int(data[i+2])
                    for j in range(nat):
                        tmol.newAtom(atype,data[i+3+j].split()[:3])
                    i+=2+nat
                else:
                    tparam[nl]+=data[i]
                i+=1
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
    return tmol,tparam

def writer(mol,f,param,coordfmt):
    """
    Save CPMD input file

    Needs both mol and param
    Respects coordfmt
    """
    for i in param:
        #always write &CPMD namelist
        if i=="&CPMD":
            f.write("&CPMD\n")
            f.write(param[i])
            f.write("&END\n")
        #expand &SYSTEM with atom and cell information
        elif i == "&SYSTEM":
            f.write("&SYSTEM\n")
            #specify coordfmt
            if coordfmt=="angstrom":
                f.write("  ANGSTROM\n")
            elif coordfmt=="crystal":
                f.write("  SCALE\n")
            elif coordfmt=="alat":
                f.write("  SCALE CARTESIAN\n")
            #cell vectors are always given explicitely
            f.write("  CELL VECTORS\n")
            vec=mol.getVec()*mol.getCellDim(coordfmt)
            for j in vec:
                f.write("  {:10.5f} {:10.5f} {:10.5f}\n".format(*j))
            #write K-Points if not gamma
            if mol.getKpoints('active')=='mpg':
                kpoints=mol.getKpoints('mpg')
                f.write("  KPOINTS MONKHORST-PACK")
                if all([float(j)==0 for j in kpoints[3:]]):
                    f.write("\n")
                else:
                    f.write("  SHIFT={:4s} {:4s} {:4s}\n".format(*kpoints[3:]))
                f.write("  {:4s} {:4s} {:4s}\n".format(*kpoints[:3]))
            elif mol.getKpoints('active')=='discrete':
                kpoints=mol.getKpoints('discrete')
                opts=mol.getKpoints('options')
                f.write("  KPOINTS")
                if opts["crystal"]:
                    f.write("  SCALED")
                if opts["bands"]:
                    f.write("  BANDS\n")
                    for j in range(0,len(kpoints),2):
                        line=kpoints[j]+kpoints[j+1]
                        f.write("  {3:4s} {0:4s} {1:4s} {2:4s} {4:4s} {5:4s} {6:4s}\n".format(*line))
                    f.write("  0 0. 0. 0. 0. 0. 0.\n")
                else:
                    f.write("\n  "+str(len(kpoints))+"\n")
                    for j in kpoints:
                        f.write("  {:4s} {:4s} {:4s} {:4s}\n".format(*j))
            #write rest of &SYSTEM namelist
            f.write(param[i])
            f.write("&END\n")
        #put coordinates and PP informations here
        elif i=="&ATOMS":
            f.write("&ATOMS\n")
            atoms=mol.getAtoms(coordfmt)
            types=mol.getTypes()
            for j in types:
                pp=mol.pse[j]['CPPP']
                if not pp:
                    pp=j+config['Default CPMD PP-suffix']
                nl=mol.pse[j]['CPNL']
                if not nl:
                    nl=config['Default CPMD Nonlocality']
                f.write("*"+pp+"\n")
                f.write("  "+nl+"\n")
                f.write("  "+str([a[0] for a in atoms].count(j))+"\n")
                for k in atoms:
                    if k[0]==j:
                        f.write("  {:10.5f} {:10.5f} {:10.5f}\n".format(*k[1]))
            f.write("&END\n")
        #write other namelists only when they're not empty
        elif i[0]=="&":
            if param[i]:
                f.write(i+"\n")
                f.write(param[i])
                f.write("&END\n")
