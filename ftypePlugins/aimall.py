name = 'AIMALL'
extension = 'sum'
argument = '-aim'

writer = None

def parser(controller,data):
    """ Parse AIMALL output to molecule

    Creates atoms from NACPs
    Other critical points will be parsed directly
    """
    controller.newMol()
    tmol = controller.getMol(-1)
    i=0
    while i<len(data):
        line = data[i].split()
        if line and line[0] == 'CP#':
            i+=1
            cptype = data[i].split()[3]
            if cptype == 'NACP':
                cptype = data[i].split()[4].rstrip('0123456789')
            tmol.create_atom(cptype,line[4:],'bohr')
        i+=1
