# -*- coding: utf-8 -*-
from collections import OrderedDict
from math import sqrt

name = 'PWScf Input'
extension = 'pwi'
argument = '-pwi'

def parser(controller,data):
    """ Parse PWScf input files

    Namelists will be parsed and saved in PW parameter set
    Supported Cards:
    - ATOMIC_SPECIES
    - ATOMIC_POSITIONS
    - K_POINTS
    - CELL_PARAMETERS
    Not supported:
    - CONSTRAINTS
    - OCCUPATIONS
    - ATOMIC_FORCES (PWSCFv5)
    """
    controller.newMol()
    tmol = controller.getMol(-1)
    controller.newPw()
    tparam = controller.getPw(-1)
    tcoord = []
    tvec = [[0,0,0],[0,0,0],[0,0,0]]
    #parse data and create tparam
    while data:
        header = data.pop(0).strip().split()
        # ignore empty lines
        if not header:
            pass
        # parse namelists
        elif header[0][0] == '&':
            tnl = OrderedDict()
            # parse entries
            line = data.pop(0).strip().split(',')
            while line[0] != '/':
                for j in range(len(line)):
                    if line[j]:
                        tnl[line[j].split('=')[0].strip()]=line[j].split('=')[1].strip()
                line = data.pop(0).strip().split(',')
            tparam[header[0].lower()]=tnl
        # parse card
        elif header[0][0].isupper():
            #ATOMIC_SPECIES:
            #Name   Weight  PP-file
            if header[0] == 'ATOMIC_SPECIES':
                for i in range(int(tparam['&system']['ntyp'])):
                    line = data.pop(0).strip().split()
                    tmol.pse[line[0]]['m'] = float(line[1])
                    tmol.pse[line[0]]['PPfile'] = line[2]
            #ATOMIC_POSITIONS fmt
            #Name   x   y   z
            elif header[0] == 'ATOMIC_POSITIONS':
                fmt = header[1].strip('{()}')
                for i in range(int(tparam['&system']['nat'])):
                    #support empty lines
                    temp = data.pop(0).strip().split()
                    while not temp:
                            temp = data.pop(0).strip().split()
                    tcoord.append(temp)
            #K_POINTS fmt
            elif header[0] == 'K_POINTS':
                #Gamma point only
                if header[1].strip('{()}') == 'gamma':
                    tmol.set_kpoints('active','gamma')
                #Monkhorst Pack Grid:
                #x y z offset
                elif header[1].strip('{()}') == 'automatic':
                    line = data.pop(0).strip().split()
                    tmol.set_kpoints('automatic',line)
                    tmol.set_kpoints('active','automatic')
                #else:
                #number of kpoints
                #x y z weight
                else:
                    nk = int(data.pop(0).strip().split()[0])
                    kpoints = []
                    for i in range(nk):
                        kpoints.append(data.pop(0).strip().split())
                    tmol.set_kpoints('disc',kpoints)
                    tmol.set_kpoints('active',header[1])
            #CELL_PARAMETERS
            #only needed if ibrav=0
            #tbd changed between pw4 and pw5, ignored for now
            elif header[0] == 'CELL_PARAMETERS':
                for i in [0,1,2]:
                    line = data.pop(0).strip().split()
                    tvec[i]=[float(x)for x in line]
            else:
                pass
    #Identify cell parameter representation, parse it
    tmol.set_celldm(tparam['&system']['celldm(1)'])
    if tparam['&system']['ibrav'] == '0':
        #check if CELL_PARAMETERS card has been read, if not present, throw error
        if tvec == [[0,0,0],[0,0,0],[0,0,0]]:
                raise ValueError('ibrav=0, but CELL_PARAMETERS missing')
        else:
                tmol.set_vec(tvec)
    elif tparam['&system']['ibrav'] == '1':
        #simple cubic
        pass
    elif tparam['&system']['ibrav'] == '2':
        #face centered cubic
        tmol.set_vec([[-0.5,0,0.5],[0,0.5,0.5],[-0.5,0.5,0]])
    elif tparam['&system']['ibrav'] == '3':
        #body centered cubic
        tmol.set_vec([[0.5,0.5,0.5],[-0.5,0.5,0.5],[-0.5,-0.5,0.5]])
    elif tparam['&system']['ibrav'] == '4':
        #hexagonal
        ca = float(tparam['&system']['celldm(3)'])
        tmol.set_vec([[1,0,0],[-0.5,sqrt(3)*0.5,0],[0,0,ca]])
    elif tparam['&system']['ibrav'] == '5':
        #trigonal
        c = float(tparam['&system']['celldm(4)'])
        tx=sqrt((1-c)/2)
        ty=sqrt((1-c)/6)
        tz=sqrt((1+2*c)/3)
        tmol.set_vec([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]])
    elif tparam['&system']['ibrav'] == '-5':
        #trigonal,alternative
        c = float(tparam['&system']['celldm(4)'])
        tx=sqrt((1-c)/2)
        ty=sqrt((1-c)/6)
        tz=sqrt((1+2*c)/3)
        u=(tz-2*sqrt(2)*ty)/sqrt(3)
        v=(tz+sqrt(2)*ty)/sqrt(3)
        tmol.set_vec([[u,v,v],[v,u,v],[v,v,u]])
    elif tparam['&system']['ibrav'] == '6':
        #simple tetragonal
        ca = float(tparam['&system']['celldm(3)'])
        tmol.set_vec([[1,0,0],[0,1,0],[0,0,ca]])
    elif tparam['&system']['ibrav'] == '7':
        #body centered tetragonal
        ca = float(tparam['&system']['celldm(3)'])
        tmol.set_vec([[0.5,-0.5,ca*0.5],[0.5,0.5,ca*0.5],[-0.5,-0.5,ca*0.5]])
    elif tparam['&system']['ibrav'] == '8':
        #simple orthorhombic
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        tmol.set_vec([[1,0,0],[0,ba,0],[0,0,ca]])
    elif tparam['&system']['ibrav'] == '9':
        #basis centered orthorhombic
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        tmol.set_vec([[0.5,ba*0.5,0],[-0.5,ba*0.5,0],[0,0,ca]])
    elif tparam['&system']['ibrav'] == '10':
        #face centered orthorhombic
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        tmol.set_vec([[0.5,0,ca*0.5],[0.5,ba*0.5,0],[0,ba*0.5,ca*0.5]])
    elif tparam['&system']['ibrav'] == '11':
        #body centered orthorhombic
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        tmol.set_vec([[0.5,ba*0.5,ca*0.5],[-0.5,ba*0.5,ca*0.5],[-0.5,-ba*0.5,ca*0.5]])
    elif tparam['&system']['ibrav'] == '12':
        #simple monoclinic
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        cg = float(tparam['&system']['celldm(4)'])
        tmol.set_vec([[1,0,0],[ba*cg,ba*sqrt(1-cg),0],[0,0,ca]])
    elif tparam['&system']['ibrav'] == '-12':
        #simple monoclinic, alternate definition
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        cb = float(tparam['&system']['celldm(5)'])
        tmol.set_vec([[1,0,0],[0,ba,0],[ca*cb,0,ca*sqrt(1-cb)]])
    elif tparam['&system']['ibrav'] == '13':
        #base centered monoclinic
        ca = float(tparam['&system']['celldm(3)'])
        ba = float(tparam['&system']['celldm(2)'])
        cg = float(tparam['&system']['celldm(4)'])
        tmol.set_vec([[0.5,0,-ca*0.5],[ba*cg,ba*sqrt(1-cg),0],[0.5,0,ca*0.5]])
    elif tparam['&system']['ibrav'] == '14':
        #base centered monoclinic
        ba = float(tparam['&system']['celldm(2)'])
        ca = float(tparam['&system']['celldm(3)'])
        cg = float(tparam['&system']['celldm(4)'])
        cb = float(tparam['&system']['celldm(5)'])
        cal = float(tparam['&system']['celldm(6)'])
        tmol.set_vec([[1,0,0],[ba*cg,ba*sqrt(1-cg),0],
                [ca*cb,ca*(cal-cb*cg)/sqrt(1-cg),ca*sqrt(1+2*cal*cb*cg-cal*cal-cb*cb-cg*cg)/sqrt(1-cg)]])
    #create atoms after creating cell:
    for i in range(len(tcoord)):
        tmol.create_atom(tcoord[i][0],tcoord[i][1:4],fmt, [int(x) for x in tcoord[i][4:]])
    #delete nat, ntype and celldm before returning to controller
    del tparam['&system']['nat']
    del tparam['&system']['ntyp']
    del tparam['&system']['celldm(1)']

def writer(mol,f,param,coordfmt):
    """
    Save PWScf input file

    Needs both mol and param
    Respects coordfmt
    """
    #&control, &system and &electron namelists are mandatory
    for i in ['&control','&system','&electrons']:
        f.write(i+'\n')
        #write all parameters
        if i == '&system':
            f.write(' nat='+str(mol.get_nat())+'\n')
            f.write(' ntyp='+str(mol.get_ntyp())+'\n')
            f.write(' celldm(1)='+str(mol.get_celldm())+'\n')
        for j in range(len(param[i])):
            f.write(' '+param[i].keys()[j]+'='+param[i].values()[j]+'\n')
        f.write('/\n\n')
    #&ions only when needed
    if '&ions' in param:
        f.write('&ions'+'\n')
        for j in range(len(param['&ions'])):
            f.write(' '+param['&ions'].keys()[j]+'='+param['&ions'].values()[j]+'\n')
        f.write('/\n\n')
    elif param['&control']['calculation'] in ["'relax'","'vc-relax'","'md'","'vc-md'"]:
        raise KeyError('&ions namelist required, but not present')
    #&cell only when needed
    if '&cell' in param:
        f.write('&cell'+'\n')
        for j in range(len(param['&cell'])):
            f.write(' '+param['&cell'].keys()[j]+'='+param['&cell'].values()[j]+'\n')
        f.write('/\n\n')
    elif param['&control']['calculation'] in ["'vc-relax'","'vc-md'"]:
        raise KeyError('&cell namelist required, but not present')
    #ATOMIC_SPECIES card:
    f.write('ATOMIC_SPECIES'+'\n')
    types = list(mol.get_types())
    for i in range(len(mol.get_types())):
        atom = types[i]
        f.write(atom+'    '+str(mol.pse[atom]['m'])+'   '+str(mol.pse[atom]['PPfile'])+'\n')
    f.write('\n')
    #ATOMIC_POSITIONS
    f.write('ATOMIC_POSITIONS'+' '+coordfmt+'\n')
    for i in range(mol.get_nat()):
        atom=mol.get_atom(i,coordfmt)
        if 0 in atom[3]:
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f} {:1d} {:1d} {:1d}'.format(
                atom[0],atom[1][0],atom[1][1],atom[1][2],atom[3][0],atom[3][1],atom[3][2])+'\n')
        else:
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
    f.write('\n')
    #K_POINTS
    f.write('K_POINTS'+' '+mol.get_kpoints('active')+'\n')
    #Gamma point only
    if mol.get_kpoints('active') == 'gamma':
        pass
    #MPGrid:
    #x y z offset
    elif mol.get_kpoints('active') == 'automatic':
        auto = mol.get_kpoints('automatic')
        f.write('{:4s} {:4s} {:4s} {:4s} {:4s} {:4s}'.format(
                auto[0],auto[1],auto[2],auto[3],auto[4],auto[5])+'\n')
    #number of kpoints
    #x y z weight
    else:
        disc=mol.get_kpoints('disc')
        f.write(str(len(disc))+'\n')
        for i in range(len(disc)):
            f.write('{:4s} {:4s} {:4s} {:4s}'.format(
                    disc[i][0],disc[i][1],disc[i][2],disc[i][3])+'\n')
    f.write('\n')
    #Cell parameters
    f.write('CELL_PARAMETERS'+'\n')
    fmt='{0[0][0]:15.10f} {0[0][1]:15.10f} {0[0][2]:15.10f}\n' + \
        '{0[1][0]:15.10f} {0[1][1]:15.10f} {0[1][2]:15.10f}\n' + \
        '{0[2][0]:15.10f} {0[2][1]:15.10f} {0[2][2]:15.10f}\n'
    f.write(fmt.format(mol.get_vec()))
