#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import numpy as np
from mol_f import set_bonds_f,make_vol_gradient

class Molecule:
    """
    Main Class for Molecule/Cell data

    Internal coordinates/sizes saved in bohr

    Includes:
    Atom symbols/coordinates (_atom_name/_atom_coord)
    Bonds between atoms (generated when requested) (_bonds)
    List of fixed atoms in PWScf calculation (_atom_fix)
    XYZ comment line (_comment)
    Cube-style volume data (_vol)
    Cell geometry (_celldm/_vec)
    K-Point settings (_kpoints)
    """

    def __init__(self,controller):
        self.pse=self._pse(controller.pse)
        self._atom_name=[]
        self._atom_coord=[]
        self._atom_fix=[]
        self._bond_cutoff=[]
        self._script_group=dict()
        self._selection=[]
        self._celldm = 1.0
        self._vec=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],'f')
        self._vecinv=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],'f')
        self._comment = ''
        self._kpoints={'active':'gamma'}

    def create_atom(self,name='C',coord=[0.,0.,0.],fmt='bohr',fix=[1,1,1]):
        """Make new atom
        name -> element symbol
        coord -> coordinates
        fmt -> str in ['bohr'/'angstrom'/'crystal'/'alat']
        fix -> 0 for no relaxation in PW
        """
        self._atom_name.append(name)
        self._atom_coord.append(self._coord_to_bohr(coord,fmt))
        self._atom_fix.append(fix)

    def insert_atom(self,pos,addat):
        """Insert copy of an atom

        pos -> id of new atom
        addat -> id of old atom
        """
        self._atom_name.insert(pos,addat[0])
        self._atom_coord.insert(pos,self._coord_to_bohr(addat[1],addat[2]))
        self._atom_fix.insert(pos,addat[3])

    def del_atom(self,index):
        """Remove atom"""
        del self._atom_name[index]
        del self._atom_coord[index]
        del self._atom_fix[index]

    ######################################################
    # SET FUNCTIONS
    ######################################################

    def set_atom(self,index,name,coord,fmt,fix=[1,1,1]):
        """Modify a given atom

        index -> atom id
        rest -> new properties
        """
        self._atom_name[index]=name
        self._atom_coord[index]=self._coord_to_bohr(coord,fmt)
        self._atom_fix[index]=fix

    def set_comment(self,comment):
        """Specify xyz comment line"""
        self._comment = comment

    def set_celldm(self,cdm,scale=False,fmt='bohr'):
        """Set new cell-dimension

        cdm -> new cell-dimension
        scale -> if True: scale coordinates of atoms (fixed crystal-coord.)
        fmt -> str in ['angstrom'/'bohr']
        """
        if fmt=='angstrom':
            cdm=cdm*1.889726125
        if scale:
            ratio=cdm/self._celldm
            for i in range(len(self._atom_coord)):
                self._atom_coord[i] = self._atom_coord[i]*ratio
        self._celldm = float(cdm)

    def set_vec(self,vec,scale=False):
        """Set new cell-vectors

        vec -> 3x3 list of new vectors
        scale -> if True: scale coordinates of atoms (fixed crystal-coord.)
        """
        vec = np.array(vec,'f')
        inv = self._vecinv
        if scale:
            for i in range(len(self._atom_coord)):
                self._atom_coord[i] = np.dot(np.dot(self._atom_coord[i],inv),vec)
        self._vec = vec
        self._vecinv = np.linalg.inv(self._vec)

    def set_kpoints(self,mode,kpoints):
        """Set and modify active k-points

        mode -> str in [active,automatic,disc]
        kpoints -> corresponding argument
        self._kpoints['active'] in [gamma,automatic,disc]
        self._kpoints['automatic'] = [x,y,z,xoff,yoff,zoff]
        self._kpoints['tpiba|crystal|tpiba_b|crystal_b'] = [[x,y,z,weight]]*nk
        """
        if mode in ['tpiba','crystal','tpiba_b','crystal_b']:
            self._kpoints['disc']=kpoints
        else:
            self._kpoints[mode]=kpoints

    #######################################################
    # COORD FMT FUNCTIONS
    # to be called only by atom set/get
    ######################################################

    def _coord_to_bohr(self,coord,fmt='bohr'):
        """Transform given coordinates to bohr

        coord -> new coordinates in given format
        fmt -> str in ['angstrom','bohr','crystal','alat']
        retval -> new coordinates in bohr
        """
        coord = np.array(coord,'f')
        if fmt == 'angstrom':
            return coord*1.889726125
        elif fmt == 'bohr':
            return coord
        elif fmt == 'crystal':
            return np.dot(coord,self._vec)*self._celldm
        elif fmt == 'alat':
            return coord*self._celldm

    def _coord_from_bohr(self,coord,fmt):
        """Transform given coordinates from bohr

        coord -> old coordinates in bohr
        fmt -> str in ['angstrom','bohr','crystal','alat']
        retval -> new coordinates in given format
        """
        if fmt == 'angstrom':
            return coord*0.52917721092
        elif fmt == 'bohr':
            return coord
        elif fmt == 'crystal':
            return np.dot(coord,self._vecinv)/self._celldm
        elif fmt == 'alat':
            return coord/self._celldm

    ######################################################
    # LOCAL PSE-DICT
    # overlay for global dict
    ######################################################

    class _pse(dict):
        def __init__(self,cpse):
            super(Molecule._pse,self).__init__()
            self.cpse=cpse

        def __getitem__(self,key):
            if not key in self:
                if key in self.cpse:
                    self[key]=self.cpse[key]
                else:
                    for i in range(len(key),0,-1):
                        if key[:i] in self.cpse:
                            self[key] = self.cpse[key[:i]]
                            break
                    if not key in self:
                        self[key]=self.cpse['X']
            return super(Molecule._pse,self).__getitem__(key)

    ######################################################
    # RETURN FUNCTIONS
    ######################################################

    def get_nat(self):
        """Return number of atoms"""
        return len(self._atom_coord)

    def get_celldm(self,fmt='bohr'):
        """Return cell-dimension

        fmt -> str in ['angstrom','bohr','crystal','alat']
        """
        if fmt=='angstrom':
            return self._celldm*0.52917721092
        else:
            return self._celldm

    def get_all_atoms(self):
        """Return names and coordinates (bohr) for all atoms"""
        return zip(self._atom_name,self._atom_coord)

    def get_atom(self,index,fmt='bohr'):
        """Return one atom

        index -> index of atom
        fmt -> str in ['angstrom','bohr','crystal','alat']
        """
        return [self._atom_name[index],self._coord_from_bohr(self._atom_coord[index],fmt),fmt,self._atom_fix[index]]

    def get_vec(self):
        """Return cell-vectors"""
        return self._vec

    def get_comment(self):
        """Return xyz-comment"""
        return self._comment

    def get_types(self):
        """Return types of atoms"""
        types = set()
        for i in self._atom_name:
            types.add(i)
        return types

    def get_ntyp(self):
        """Return number of types of atoms"""
        return len(self.get_types())

    def get_center(self):
        """Return center-coordinates of cell"""
        return (self._vec[0]+self._vec[1]+self._vec[2])*self._celldm/2

    def get_kpoints(self,mode):
        """Return active k-points or k-point-settings"""
        if mode in ['crystal','tpiba','crystal_b','tpiba_b']:
            return self._kpoints['disc']
        else:
            return self._kpoints[mode]

    ######################################################
    # SELECTION FUNCTIONS
    ######################################################

    def add_selection(self,idx,offset):
        """Add an item to selection

        idx -> index of atom
        offset -> offset in multiples of vec
        """
        item = [idx,offset]
        if item in self._selection:
            self._selection.remove(item)
        else:
            self._selection.append(item)

    def del_selection(self,idx=None,offset=None):
        """Delete one or all atoms from selection

        If both idx and offset are given, specific
        periodic image of atom will be deleted,
        else the whole selection will be cleared.
        """
        if idx and offset:
            self._selection.remove([idx,offset])
        elif not idx and not offset:
            self._selection = []

    def get_selection(self):
        """Return the selection"""
        return self._selection

    ######################################################
    # BOND FUNCTIONS
    ######################################################

    def get_bonds(self):
        """Return bonds

        Sets bonds if not present.
        Returns list of lists.
        First entry: Bonds inside of cell
        Other entries: Periodic bonds (see set_bonds)
        """
        if not hasattr(self,'_bonds'): self.set_bonds()
        return self._bonds

    def set_bonds(self):
        """Set bonds

        Cutoff criteria:
        Maximum bondlength: 1.85 Å
        Bondlength for hydrogen: 1.2 Å
        Generates list of lists, each lists contains:
        [idx_at1,idx_at2,offset_at1,offset_at2,distance]
        Actual setting performed by fortran-routine set_bonds_f
        """
        self._bonds=[[],[],[],[],[],[],[],[]]
        if len(self._atom_coord)<2:
            return
        at_c = self._atom_coord
        cutoff=np.array([3.5]*len(at_c),'f')
        for i in range(len(at_c)):
            if self._atom_name[i] == 'H':
                cutoff[i]=2.27
        n=np.zeros(3)
        v = self.get_vec()*self.get_celldm()
        off = [[(n,n)],                             #orig
                [(v[0],n)],                         #x
                [(v[1],n)],                         #y
                [(v[0]+v[1],n),(v[0],v[1])],        #xy,x-y
                [(v[2],n)],                         #z
                [(v[0]+v[2],n),(v[0],v[2])],        #xz,x-z
                [(v[1]+v[2],n),(v[1],v[2])],        #yz,y-z
                [(v[0]+v[1]+v[2],n),(v[0]+v[1],v[2]),(v[0]+v[2],v[1]),(v[1]+v[2],v[0])]] #xyz,xy-z,x-yz,-xyz
        for k,os in enumerate(off):
            for i in os:
                nbnds,at1,at2,dist = set_bonds_f(at_c,cutoff,i)
                self._bonds[k].extend(zip(at1[:nbnds],at2[:nbnds],[i]*nbnds,dist[:nbnds]))

    #####################################################
    # EDIT FUNCTIONS
    #####################################################

    def mult(self,x,y,z):
        """Multiply cell

        x,y,z -> Integer multipliers along corresponding vector
        """
        vec = self._vec*self._celldm
        mult = [x,y,z]
        for k in range(3):
            nat = self.get_nat()
            self._atom_name = mult[k]*self._atom_name
            self._atom_fix = mult[k]*self._atom_fix
            self._atom_coord = mult[k]*self._atom_coord
            for i in range(1,mult[k]):
                for j in range(i*nat,(i+1)*nat):
                    self._atom_coord[j]=self._atom_coord[j]+i*vec[k]
        self.set_vec(self._vec*[[x],[y],[z]])

    def crop(self):
        """Crop all atoms outside of cell"""
        nat=self.get_nat()
        dellist = []
        for i in range(nat):
            at=self.get_atom(i,'crystal')
            if np.any(at[1]>1.000001) or np.any(at[1]<-0.000001) or np.any(np.isclose(at[1]-1,0,atol=1.e-6)):
                dellist.append(i)
        dellist.reverse()
        for i in dellist:
            self.del_atom(i)

    def wrap(self):
        """Wrap all atoms outside of cell to the inside"""
        nat = self.get_nat()
        for i in range(nat):
            at=self.get_atom(i,'crystal')
            self.set_atom(i,at[0],at[1]%1,'crystal')

    def reshape(self,newvec):
        """Reshape cell

        newvec -> new cell vectors

        Transforms cell so new cell is a cut of
        bulk-material with given geometry.
        PBC must be manually conserved!
        """
        self.wrap()
        newdim=abs(np.array(newvec)).sum(0)
        olddim=abs(self._vec).sum(0)
        m = 1
        while np.any(olddim*m<newdim):
            m+=1
        self.mult(m,m,m)
        oldcenter=self.get_center()
        self.set_vec(newvec)
        newcenter=self.get_center()
        for i in range(self.get_nat()):
            self._atom_coord[i]+=newcenter-oldcenter
        self.crop()

    def align(self,vec,direc):
        """Align cell vectors

        vec -> int specifying the vector to align
        direc -> str specifying the direction to align to
        """

        if direc == 'x':
            d = np.array([1,0,0],'f')
        elif direc == 'y':
            d = np.array([0,1,0],'f')
        elif direc == 'z':
            d = np.array([0,0,1],'f')
        else:
            raise ValueError('Align vectors: invalid direction!')

        v = self.get_vec()[int(vec)]
        v = v/np.linalg.norm(v)

        if np.all(np.equal(abs(v),d)):
            return

        ax = np.cross(v,d)
        ax = ax/np.linalg.norm(ax)
        c = np.float(np.dot(v,d))
        theta = np.arccos(np.dot(v,d))
        s=np.float(-np.sin(theta))
        ic=np.float(1.-c)
        mat=np.array([[ic*ax[0]*ax[0]+c,ic*ax[0]*ax[1]-s*ax[2],ic*ax[0]*ax[2]+s*ax[1]],
                      [ic*ax[0]*ax[1]+s*ax[2],ic*ax[1]*ax[1]+c,ic*ax[1]*ax[2]-s*ax[0]],
                      [ic*ax[0]*ax[2]-s*ax[1],ic*ax[1]*ax[2]+s*ax[0],ic*ax[2]*ax[2]+c]],'f')
        mat = np.array([np.dot(i,mat) for i in self.get_vec()])
        self.set_vec(mat,scale=True)

    #####################################################
    # VOLUME DATA FUNCTIONS
    #####################################################

    def set_vol(self,dim,vol,off):
        """Set volume data

        dim -> list of dimension of data-grid
        vol -> data-grid
        off -> offset of data

        Parses string-list containing cube-style volume data
        to list of shape dim[0]*dim[1]*dim[2]
        """
        self._vol=np.zeros(dim,'f')
        self._vol_off = np.array(off,'f')
        i=0
        j=0
        line=vol[i].split()
        for x in range(dim[0]):
            for y in range(dim[1]):
                for z in range(dim[2]):
                    self._vol[x][y][z]=float(line[j])
                    j+=1
                    if j==len(line) and i<(len(vol)-1):
                        j=0
                        i+=1
                        line=vol[i].split()

    def get_vol(self):
        """Return volume data"""
        return self._vol

    def get_vol_offset(self):
        """Return offset of volume data"""
        return self._vol_off

    def set_vol_gradient(self):
        """Calculate the gradient of volume data"""
        self._vol_grad = make_vol_gradient(self._vol)
        #grad_x = np.roll(self._vol,-1,0) - np.roll(self._vol,1,0)
        #grad_y = np.roll(self._vol,-1,1) - np.roll(self._vol,1,1)
        #grad_z = np.roll(self._vol,-1,2) - np.roll(self._vol,1,2)
        #self._vol_grad = zip(grad_x,grad_y,grad_z)

    def get_vol_gradient(self):
        """Return volume gradient"""
        if not hasattr(self,'_vol_grad'):
            self.set_vol_gradient()
        return self._vol_grad

    #####################################################
    # SCRIPTING SUPPORT
    #####################################################

    def evalScript(self,script):
        """Evaluate script

        script -> list of str containing script

        Operators are defined in ops-dictionary.
        Defined groups of atoms are saved in instanced list '_script_group'.
        Other evaluated operations are placed on stack and
        executed after successful parsing.
        """
        script=script.strip().replace('\n',' ')
        rep=1
        #dictionary returns closure containing target operation and argument-check
        ops={'rot':self._evalArgs(self._rotate,'lavo'),
                'shi':self._evalArgs(self._shift,'lv'),
                'def':self._evalArgs('define','ls'),
                'mir':self._evalArgs(self._mirror,'lvvo'),
                'rep':self._evalArgs('repeat','i')}
        stack=[]
        #check for errors, parse and prepare
        while script:
            try:
                op,script = script.split(' ',1)
                op,script = ops[op[0:3]](script)
            except KeyError as e:
                return 'Wrong Op: '+script[0]
            except StandardError as e:
                print(e.message)
                return e.message
            else:
                if op[0]=='define':
                    self._script_group[op[2]]=op[1]
                elif op[0]=='repeat':
                    rep=op[1]
                else:
                    for i in range(rep):
                        stack.append(op)
                    rep=1
        # delete previous definitions
        self._script_group={}
        #if everything went well, execute operations
        for op in stack:
            op[0](*op[1:])
        return 'Success!'

    def _evalArgs(self,op,args):
        """Evaluate script-op-arguments

        op -> Given operator
        args -> required arguments for operator
        retval -> closure for evaluation of given arguments
        """
        def getArg(arglist,sep=' '):
            arglist = arglist.split(sep,1)
            arg = arglist[0]
            if sep!=' ':
                arg+=sep
            if len(arglist)>1:
                rest = arglist[1].strip()
            else:
                rest = ''
            return arg,rest

        def evArgs(arglist):
            res = [op]
            for t in args:
                #list of atoms
                if t == 'l':
                    #list of indices
                    if arglist[0] == '[':
                        arg,arglist=getArg(arglist,']')
                        arg = arg.strip('[]').split(',')
                        l=[]
                        for j in arg:
                            #interpret range
                            if '-' in j:
                                low,high = j.split('-')
                                l.extend(range(int(low),int(high)+1))
                            #atom index
                            else:
                                l.append(int(j))
                        if not np.all(np.less(l,len(self._atom_coord))):
                            raise IndexError('Atom index out of range')
                        res.append(l)
                    else:
                        arg,arglist = getArg(arglist)
                        #reference to definition
                        if arg.isalpha():
                            if arg == 'all':
                                res.append(range(len(self._atom_name)))
                            elif arg == 'sel':
                                res.append(set([a[0] for a in self._selection]))
                            else:
                                res.append(self._script_group[arg])
                        #single index
                        elif arg.isdigit():
                            res.append([int(arg)])
                        elif '-' in arg:
                            low,high = arg.split('-')
                            res.append(range(int(low),int(high)+1))
                #angle
                elif t == 'a':
                    arg,arglist=getArg(arglist)
                    res.append(float(arg))
                #index for loops
                elif t == 'i':
                    arg,arglist=getArg(arglist)
                    res.append(int(arg))
                #arbitrary names for defined groups
                elif t == 's':
                    arg,arglist=getArg(arglist)
                    res.append(arg)
                #valid vector for operations
                elif t in 'vo':
                    #nothing left and optional:
                    if not arglist and t == 'o':
                        continue
                    #explicit vector
                    if arglist[0] == '(':
                        arg,arglist = getArg(arglist,')')
                        arg = eval(arg)
                        if len(arg)==4 and type(arg[3]) is str:
                            arg = self._coord_to_bohr(arg[0:3],arg[3])
                        elif len(arg)==3:
                            arg = self._coord_to_bohr(arg)
                        else:
                            raise ValueError('Not a valid vector: '+str(arg))
                    #implicit vector
                    elif arglist[0].isdigit():
                        arg,arglist = getArg(arglist)
                        #difference between atoms
                        if '-' in arg:
                            arg=arg.split('-')
                            arg=self._atom_coord[int(arg[0])]-self._atom_coord[int(arg[1])]
                        #position of atom
                        else:
                            arg=np.array(self._atom_coord[int(arg)])
                    #fail when not a vector and vector needed
                    elif t=='v':
                        raise ValueError('Not a valid vector: '+str(arg))
                    #continue when not a vector and not needed
                    else:
                        continue
                    res.append(arg)
            return res,arglist
        return evArgs

    def _rotate(self,atoms,angle,ax,shift=np.zeros(3)):
        """Rotate group of atoms

        atoms -> list of atoms
        angle -> angle in degree
        ax -> rotation-axis
        shift -> optional vector shifting the rotation-axis
        """
        angle=np.radians(angle)
        c=np.float(np.cos(angle))
        s=np.float(-np.sin(angle))
        ic=np.float(1.-c)
        ax=ax/np.linalg.norm(ax)
        mat=np.array([[ic*ax[0]*ax[0]+c,ic*ax[0]*ax[1]-s*ax[2],ic*ax[0]*ax[2]+s*ax[1]],
                      [ic*ax[0]*ax[1]+s*ax[2],ic*ax[1]*ax[1]+c,ic*ax[1]*ax[2]-s*ax[0]],
                      [ic*ax[0]*ax[2]-s*ax[1],ic*ax[1]*ax[2]+s*ax[0],ic*ax[2]*ax[2]+c]],'f')
        for i in atoms:
            self._atom_coord[i]=np.dot(self._atom_coord[i]-shift,mat)+shift

    def _shift(self,atoms,vector):
        """Shift group of atoms

        atoms -> list of atoms
        vector -> shift-vector
        """
        for i in atoms:
            self._atom_coord[i]+=np.array(vector,'f')

    def _mirror(self,atoms,v1,v2,shift=np.zeros(3)):
        """Mirror group of atoms

        atoms -> group of atoms
        v1,v2 -> vectors defining the plane
        shift -> optional vector shifting the mirror-plane
        """
        normal=np.cross(v1,v2)
        normal=normal/np.linalg.norm(normal)
        for i in atoms:
            pos=self._atom_coord[i]
            proj = np.dot(pos-shift,normal)*normal
            self._atom_coord[i]=pos-2*proj

