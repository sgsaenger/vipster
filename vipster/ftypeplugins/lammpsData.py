# -*- coding: utf-8 -*-
from ..molecule import Molecule
from collections import OrderedDict

name = "Lammps Data"
extension = 'lmp'
argument = 'lmp'
param={"default":{"type":"lmp",
                  "name":"New LAMMPS",
                  "atom_style":"atomic",
                  "bonds":False,"angles":False,
                  "dihedrals":False,"impropers":False}}
lammps_atom_style=OrderedDict([
                ("angle",     "{0:d} {1:d} {2:d} {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("atomic",    "{0:d} {2:d} {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("body",      "{0:d} {2:d} 0 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("bond",      "{0:d} {1:d} {2:d} {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("charge",    "{0:d} {2:d} 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("dipole",    "{0:d} {2:d} 0 {3:15.10f} {4:15.10f} {5:15.10f} 0 0 0\n"),
                ("electron",  "{0:d} {2:d} 0 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("ellipsoid", "{0:d} {2:d} 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("full",      "{0:d} {1:d} {2:d} 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("line",      "{0:d} {1:d} {2:d} 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("meso",      "{0:d} {2:d} 0 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("molecular", "{0:d} {1:d} {2:d} {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("peri",      "{0:d} {2:d} 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("smd",       "{0:d} {2:d} {1:d} 0 0 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("sphere",    "{0:d} {2:d} 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("template",  "{0:d} {1:d} 0 0 {2:d} {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("tri",       "{0:d} {1:d} {2:d} 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n"),
                ("wavepacket","{0:d} {2:d} 0 0 0 0 0 0 {3:15.10f} {4:15.10f} {5:15.10f}\n")])

def parser(name,data):
    """ Parse Lammps data file

    Element parsed via comment in 'Masses' part
    Only orthogonal cells supported
    Assumes angstrom and atom_style != dipole/hybrid
    """
    tmol = Molecule(name)
    i=0
    tvec=[[0,0,0],[0,0,0],[0,0,0]]
    while i< len(data):
        line = data[i].strip()
        if 'atoms' in line:
            nat = int(line.split()[0])
        elif 'atom types' in line:
            types = [0]*int(line.split()[0])
        elif 'xlo xhi' in line:
            tvec[0][0] = float(line.split()[1])-float(line.split()[0])
        elif 'ylo yhi' in line:
            tvec[1][1] = float(line.split()[1])-float(line.split()[0])
        elif 'zlo zhi' in line:
            tvec[2][2] = float(line.split()[1])-float(line.split()[0])
            tmol.setVec(tvec)
            tmol.setCellDim(1,fmt='angstrom')
        elif 'Masses' in line:
            for j in range(i+2,i+2+len(types)):
                if '#' in data[j]:
                    types[int(data[j].split()[0])-1]=data[j].split('#')[1].strip()
                else:
                    raise NotImplementedError('cannot assign elements via masses')
            i+=len(types)+1
        elif 'Atoms' in line:
            #guess where atom-type is found:
            atomlist=[]
            for j in range(i+2,i+2+nat):
                atomlist.append(data[j].strip().split())
            nargs=len(atomlist[0])
            atypepos=1
            for j in range(1,nargs):
                try:
                    int(atomlist[0][j])
                except:
                    continue
                atomids = sorted(set([int(k[j]) for k in atomlist]))
                if atomids==range(1,len(types)+1):
                    atypepos=j
                    break
            tmol.newAtoms(nat)
            for j in range(nat):
                at = atomlist[j]
                tmol.setAtom(j,types[int(at[atypepos])-1],at[-3:])
            tmol.setFmt('angstrom',scale=True)
        i+=1
    return tmol,None

def writer(mol,f,param):
    """
    Save Lammps input data

    Cell parameters have to be in lower triangular form
    Non orthogonal cell not supported for now
    Assumes angstrom
    """
    #calculate necessary values:
    if param["bonds"] or param["angles"]\
            or param["dihedrals"] or param["impropers"]\
            or param["atom_style"] in ["angle","bond","full","line","molecular","smd","template","tri"]:
        bondlist=[]
        for i in mol.getBonds():
            for j in i:
                bondlist.append(tuple(sorted(j[0:2])))
        bondlist=sorted(set(bondlist))
        bondtypes=[]
        for i in bondlist:
            bondtypes.append(tuple(sorted([mol.getAtom(j)[0] for j in i])))
        bondtypelist=list(set(bondtypes))
    if param["angles"] or param["dihedrals"] or param["impropers"]:
        anglelist=[]
        for i in range(len(bondlist)):
            for j in range(i+1,len(bondlist)):
                if bondlist[i][0]==bondlist[j][0]:
                    anglelist.append(tuple([bondlist[i][1],bondlist[i][0],bondlist[j][1]]))
                elif bondlist[i][1]==bondlist[j][0]:
                    anglelist.append(tuple([bondlist[i][0],bondlist[i][1],bondlist[j][1]]))
                elif bondlist[i][1]==bondlist[j][1]:
                    anglelist.append(tuple([bondlist[i][0],bondlist[i][1],bondlist[j][0]]))
        angletypes=[]
        for i in anglelist:
            at=[mol.getAtom(j)[0] for j in i]
            if at[0]>at[2]:
                at.reverse()
            angletypes.append(tuple(at))
        angletypelist=list(set(angletypes))
    if param["dihedrals"]:
        dihedrallist=[]
        for i in range(len(anglelist)):
            for j in range(i+1,len(anglelist)):
                a1=anglelist[i]
                a2=anglelist[j]
                if a1[1]==a2[0] and a2[1]==a1[2]:
                    dihedrallist.append(tuple([a1[0],a1[1],a1[2],a2[2]]))
                elif a1[1]==a2[0] and a2[1]==a1[0]:
                    dihedrallist.append(tuple([a1[2],a1[1],a1[0],a2[2]]))
                elif a1[1]==a2[2] and a2[1]==a1[2]:
                    dihedrallist.append(tuple([a1[0],a1[1],a1[2],a2[0]]))
                elif a1[1]==a2[2] and a2[1]==a1[0]:
                    dihedrallist.append(tuple([a1[2],a1[1],a1[0],a2[0]]))
        dihedraltypes=[]
        for i in dihedrallist:
            dt=[mol.getAtom(j)[0] for j in i]
            if dt[0]>dt[3]:
                dt.reverse()
            elif dt[0]==dt[3] and dt[1]>dt[2]:
                dt.reverse()
            dihedraltypes.append(tuple(dt))
        dihedraltypelist=list(set(dihedraltypes))
    if param["impropers"]:
        from itertools import combinations
        improperlist=[]
        for i in range(mol.nat):
            neigh=[]
            for j in bondlist:
                if i in j:
                    neigh.append(j[not j.index(i)])
            for j in combinations(neigh,3):
                improperlist.append((i,)+j)
        impropertypes=[]
        for i in improperlist:
            it=[mol.getAtom(j)[0] for j in i]
            if it[0]>it[3]:
                it.reverse()
            elif it[0]==dt[3] and dt[1]>dt[2]:
                it.reverse()
            impropertypes.append(tuple(dt))
        impropertypelist=list(set(impropertypes))
    moleculeid = [0]*mol.nat
    if param["atom_style"] in ["angle","bond","full","line","molecular","smd","template","tri"]:
        molbondlist = [set(i) for i in bondlist]
        def groupsets(setlist):
            tlist=[setlist[0]]
            for i in setlist[1:]:
                matched=False
                for j in tlist:
                    if i&j:
                        j.update(i)
                        matched=True
                if not matched:
                    tlist.append(i)
            return tlist
        moleculelist = groupsets(molbondlist)
        for i,m in enumerate(moleculelist):
            for at in m:
                moleculeid[at]=i
    #write header:
    f.write('\n'+str(mol.nat)+' atoms\n')
    f.write(str(mol.ntyp)+' atom types\n')
    if param["bonds"] and bondlist:
        f.write(str(len(bondlist))+' bonds\n')
        f.write(str(len(bondtypelist))+' bond types\n')
        for i,j in enumerate(bondtypelist):
            f.write('#{:d} {:s} {:s}\n'.format(i+1,*j))
    if param["angles"] and anglelist:
        f.write(str(len(anglelist))+' angles\n')
        f.write(str(len(angletypelist))+' angle types\n')
        for i,j in enumerate(angletypelist):
            f.write('#{:d} {:s} {:s} {:s}\n'.format(i+1,*j))
    if param["dihedrals"] and dihedrallist:
        f.write(str(len(dihedrallist))+' dihedrals\n')
        f.write(str(len(dihedraltypelist))+' dihedral types\n')
        for i,j in enumerate(dihedraltypelist):
            f.write('#{:d} {:s} {:s} {:s} {:s}\n'.format(i+1,*j))
    if param["impropers"] and improperlist:
        f.write(str(len(improperlist))+' impropers\n')
        f.write(str(len(impropertypelist))+' improper types\n')
        for i,j in enumerate(impropertypelist):
            f.write('#{:d} {:s} {:s} {:s} {:s}\n'.format(i+1,*j))
    f.write('\n')
    #if cell is orthogonal, write vectors:
    v=mol.getVec()*mol.getCellDim(fmt='angstrom')
    if not v.diagonal(1).any() and not v.diagonal(2).any():
        if not v.diagonal(-1).any() and not v.diagonal(-2).any():
            f.write('{:.5f} {:.5f} xlo xhi\n'.format(0.0,v[0][0]))
            f.write('{:.5f} {:.5f} ylo yhi\n'.format(0.0,v[1][1]))
            f.write('{:.5f} {:.5f} zlo zhi\n\n'.format(0.0,v[2][2]))
        else:
            raise TypeError('Non-orthogonal cell not yet supported')
    else:
        raise TypeError('Cell not in suitable Lammps format')
    #Masses section: (always)
    f.write('Masses\n\n')
    atomtypes=list(mol.getTypes())
    for i,j in enumerate(atomtypes):
        f.write('{:d} {:2.4f} #{:s}\n'.format(i+1,mol.pse[j]['m'],j))
    #Atoms section: (always)
    f.write('\nAtoms\n\n')
    string=lammps_atom_style[param["atom_style"]]
    for i in range(mol.nat):
        at=mol.getAtom(i,fmt='angstrom')
        f.write(string.format(i+1,moleculeid[i]+1,atomtypes.index(at[0])+1,*at[1]))
    #Bonds section:
    convertIndices=lambda x: list(map(lambda y: tuple(map(lambda z: z+1,y)),x))
    if param["bonds"] and bondlist:
        f.write('\nBonds\n\n')
        bondlist=convertIndices(bondlist)
        for i,j in enumerate(bondlist):
            f.write('{:d} {:d} {:d} {:d}\n'.format(i+1,bondtypelist.index(bondtypes[i])+1,*j))
    #Angles section:
    if param["angles"] and anglelist:
        f.write('\nAngles\n\n')
        anglelist=convertIndices(anglelist)
        for i,j in enumerate(anglelist):
            f.write('{:d} {:d} {:d} {:d} {:d}\n'.format(i+1,angletypelist.index(angletypes[i])+1,*j))
    #Dihedrals section:
    if param["dihedrals"] and dihedrallist:
        f.write('\nDihedrals\n\n')
        dihedrallist=convertIndices(dihedrallist)
        for i,j in enumerate(dihedrallist):
            f.write('{:d} {:d} {:d} {:d} {:d} {:d}\n'.format(i+1,dihedraltypelist.index(dihedraltypes[i])+1,*j))
    #Impropers section:
    if param["impropers"] and improperlist:
        f.write('\nImpropers\n\n')
        improperlist=convertIndices(improperlist)
        for i,j in enumerate(improperlist):
            f.write('{:d} {:d} {:d} {:d} {:d} {:d}\n'.format(i+1,impropertypelist.index(impropertypes[i])+1,*j))
    f.write('\n')
