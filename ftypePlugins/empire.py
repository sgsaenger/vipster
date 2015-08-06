name = 'Empire XYZ'
extension = 'xyz'
argument = '-emp'

def parser(controller,data):
    """ Parse Empire specific xyz file """
    pass

def writer(mol,f,param,coordfmt):
    """
    Save Empire input file

    Same as xyz, cell parameters given after coordinates
    """
    # fixed format nat and comment
    f.write(str(mol.get_nat())+'\n')
    f.write('Hamil=PM3 calc=spt Periodic\n')
    # write coordinates
    for cntat in range(0,mol.get_nat()):
            atom=mol.get_atom(cntat,'angstrom')
            f.write('{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                         atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
    f.write('\n')
    vec = mol.get_vec()*mol.get_celldm()*0.52917721092
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[0]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[1]))
    f.write('{:.10f} {:.10f} {:.10f}\n'.format(*vec[2]))
