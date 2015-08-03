name = 'xyz'
extension = 'xyz'
argument = '-xyz'

writer = None

def parser(controller,data):
    """ Parse xyz files (angstrom) """
    # create list of mol, trajectory support
    controller.newTrajectory()
    tmol = controller.getMol(-1)
    i=0
    while i < len(data):
            # handle empty lines at eof or between molecules
            if not data[i].strip().isdigit():
                    i+=1
                    continue
            #fixed format nat and comment
            tmol.newMol()
            nat = int(data[i])
            tmol._comment = data[i+1].strip()
            #read coordinates and types
            for j in range(i+2,i+nat+2):
                    line = data[j].split()
                    tmol.create_atom(line[0],map(float,line[1:4]),'angstrom')
            i+=nat+2
