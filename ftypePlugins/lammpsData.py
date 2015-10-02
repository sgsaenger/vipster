name = "Lammps Data file"
extension = 'lmp'
argument = '-lmp'

def parser(controller,data):
    """ Parse Lammps data file

    Preliminary implementation!
    Element parsed via comment in 'Masses' part
    Only orthogonal cells supported
    Assumes angstrom
    """
    controller.newMol()
    tmol=controller.getMol(-1)
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
            tmol.set_vec(tvec)
            tmol.set_celldm(1,fmt='angstrom')
        elif 'Masses' in line:
            for j in range(i+2,i+2+len(types)):
                if '#' in data[j]:
                    types[int(data[j].split()[0])-1]=data[j].split('#')[1].strip()
                else:
                    raise NotImplementedError('cannot assign elements via masses')
            i+=len(types)+1
        elif 'Atoms' in line:
            for j in range(i+2,i+2+nat):
                at = data[j].strip().split()
                tmol.create_atom(types[int(at[1])-1],at[-3:],'angstrom')
        i+=1

def writer(mol,f,param,coordfmt):
    """
    Save Lammps input data

    Preliminary implementation!
    Cell parameters have to be in lower triangular form
    Non orthogonal cell not supported for now
    Assumes angstrom until dialog is established
    """
    f.write('\n'+str(mol.get_nat())+' atoms\n')
    f.write(str(mol.get_ntyp())+' atom types\n\n')
    #check if box is orthogonal:
    v=mol.get_vec()*mol.get_celldm(fmt='angstrom')
    if not v.diagonal(1).any() and not v.diagonal(2).any():
        if not v.diagonal(-1).any() and not v.diagonal(-2).any():
            f.write('{:.5f} {:.5f} xlo xhi\n'.format(0.0,v[0][0]))
            f.write('{:.5f} {:.5f} ylo yhi\n'.format(0.0,v[1][1]))
            f.write('{:.5f} {:.5f} zlo zhi\n\n'.format(0.0,v[2][2]))
        else:
            raise TypeError('Non-orthogonal cell not yet supported')
    else:
        raise TypeError('Cell not in suitable Lammps format')
    f.write('Masses\n\n')
    t=list(mol.get_types())
    for i,j in enumerate(t):
        f.write('{:d} {:2.4f} #{:s}\n'.format(i+1,mol.pse[j]['m'],j))
    f.write('\nAtoms\n\n')
    for i in range(mol.get_nat()):
        at=mol.get_atom(i,'angstrom')
        f.write(('{:d} {:d}'+' {:d}'*1+' {:15.10f} {:15.10f} {:15.10f}\n').format(
            i,t.index(at[0])+1,0,*at[1]))
    f.write('\n')
