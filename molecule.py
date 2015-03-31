#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

from time import time
import numpy as np
from mol_f import set_bonds_f

######################################################################
# MOLECULE CLASS
######################################################################
class Molecule:

        def __init__(self):
                # set atom list
                self._atom_name=[]
                self._atom_coord=[]
                self._script_group=dict()
                self._celldm = 1.0
                self._vec=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self._vecinv=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                self._comment = ''

        ######################################################
        # ATOM FUNCTIONS
        ######################################################

        # append new atom
        def create_atom(self,name='C',x=0.,y=0.,z=0.,fmt='angstrom'):
                self._atom_name.append(name)
                self._atom_coord.append(self._set_coord([x,y,z],fmt))

        # append copy of existing atom
        def append_atom_cp(self,addat):
                self._atom_name.append(addat[0])
                self._atom_coord.append(self._set_coord(addat[1],addat[2]))

        # insert atom at given position
        def insert_atom(self,pos,addat):
                self._atom_name.insert(pos,addat[0])
                self._atom_coord.insert(pos,self._set_coord(addat[1],addat[2]))

        # remove atom
        def del_atom(self,index):
                del self._atom_name[index]
                del self._atom_coord[index]

        ######################################################
        # SET FUNCTIONS
        ######################################################

        def set_atom(self,index,name,coord,fmt):
                self._atom_name[index]=name
                self._atom_coord[index]=self._set_coord(coord,fmt)

        def set_comment(self,comment):
                self._comment = comment

        # set celldm
        def set_celldm(self,cdm):
                self._celldm = float(cdm)

        # set vectors in bohr
        def set_vec(self,vec):
                self._vec = np.array(vec)
                self._vecinv = np.linalg.inv(self._vec)

        #######################################################
        # COORD FMT FUNCTIONS
        # to be called only by atom set/get
        ######################################################

        def _set_coord(self,coord,fmt='bohr'):
                coord = np.array(coord,'f')
                if fmt == 'angstrom':
                        return coord/0.52917721092
                elif fmt == 'bohr':
                        return coord
                elif fmt == 'crystal':
                        return np.dot(coord,self._vec)*self._celldm
                elif fmt == 'alat':
                        return coord*self._celldm

        def _get_coord(self,coord,fmt):
                if fmt == 'angstrom':
                        return coord*0.52917721092
                elif fmt == 'bohr':
                        return coord
                elif fmt == 'crystal':
                        return np.dot(coord,self._vecinv)/self._celldm
                elif fmt == 'alat':
                        return coord/self._celldm

        ######################################################
        # RETURN FUNCTIONS
        ######################################################
        def get_nat(self):
                return len(self._atom_coord)

        def get_celldm(self):
                return self._celldm

        def get_all_atoms(self):
                return zip(self._atom_name,self._atom_coord)

        def get_atom(self,index,fmt='bohr'):
                return [self._atom_name[index],self._get_coord(self._atom_coord[index],fmt),fmt]

        def get_vec(self):
                return self._vec

        def get_comment(self):
                return self._comment

        def get_types(self):
                types = set()
                for i in self._atom_name:
                    types.add(i)
                return types

        def get_ntyp(self):
                return len(self.get_types())

        def get_center(self):
                return (self._vec[0]+self._vec[1]+self._vec[2])*self._celldm/2

        ######################################################
        # BOND FUNCTIONS
        ######################################################
        #TODO: unify _bonds and _pbc_bonds, just return _pbc_bonds[0]

        def get_bonds(self):
                if not hasattr(self,'_bonds'): self.set_bonds()
                return self._bonds

        def get_pbc_bonds(self):
                if not hasattr(self,'_pbc_bonds'): self.set_pbc_bonds()
                return self._pbc_bonds

        def set_bonds(self):
                self.set_bonds_2(0)

        def set_bonds_2(self,pbc):
                if not pbc:
                        if len(self._atom_coord)>1:
                                at_c = self._atom_coord
                                at_n = self._atom_name
                                n = np.zeros(3)
                                print (n,n)
                                nbnds,b1,b2,n1,n2 = set_bonds_f(at_c,at_n,(n,n))
                                self._bonds=zip(b1[:nbnds],b2[:nbnds],n1[:nbnds],n2[:nbnds])
                        else:
                                self._bonds = []
                else:
                        if len(self._atom_coord)>1:
                                self._pbc_bonds=[[],[],[],[],[],[],[],[]]
                                at_c = self._atom_coord
                                at_n = self._atom_name
                                v = self.get_vec()*self.get_celldm()
                                n = np.zeros(3)
                                off = [[(n,n)],                             #orig
                                        [(v[0],n)],                         #x
                                        [(v[1],n)],                         #y
                                        [(v[0]+v[1],n),(v[0],v[1])],        #xy,x-y
                                        [(v[2],n)],                         #z
                                        [(v[0]+v[2],n),(v[0],v[2])],        #xz,x-z
                                        [(v[1]+v[2],n),(v[1],v[2])],        #yz,y-z
                                        [(v[0]+v[1]+v[2],n),(v[0]+v[1],v[2]),(v[0]+v[2],v[1]),(v[1]+v[2],v[0])]] #xyz,xy-z,x-yz,-xyz
                                for k,os in enumerate(off):
                                        self._pbc_bonds[k]=[]
                                        for i in os:
                                                print i
                                                nbonds,b1,b2,n1,n2 = set_bonds_f(at_c,at_n,i)
                                                self._pbc_bonds[k].extend(zip(b1[:nbonds],b2[:nbonds],n1[:nbonds],n2[:nbonds]))
                        else:
                                self._pbc_bonds=[self._bonds,[],[],[],[],[],[],[]]
                print self._bonds
                print self._pbc_bonds

        def set_pbc_bonds(self):
                self.set_bonds_2(1)
                return
                nat = self.get_nat()
                self._pbc_bonds=[self.get_bonds(),[],[],[],[],[],[],[]]
                v = self.get_vec()*self.get_celldm()
                off = [[(0,0)],                             #orig
                        [(v[0],0)],                         #x
                        [(v[1],0)],                         #y
                        [(v[0]+v[1],0),(v[0],v[1])],        #xy,x-y
                        [(v[2],0)],                         #z
                        [(v[0]+v[2],0),(v[0],v[2])],        #xz,x-z
                        [(v[1]+v[2],0),(v[1],v[2])],        #yz,y-z
                        [(v[0]+v[1]+v[2],0),(v[0]+v[1],v[2]),(v[0]+v[2],v[1]),(v[1]+v[2],v[0])]] #xyz,xy-z,x-yz,-xyz
                at_c = self._atom_coord
                at_n = self._atom_name
                for i in range(nat):
                        for j in range(nat):
                                dist_at = self._atom_coord[i] - self._atom_coord[j]
                                for k in [1,2,3,4,5,6,7]:
                                    for l in off[k]:
                                        dist= dist_at+l[0]-l[1]

                                        #cancel if distance in one direction is greater than allowed bond length
                                        if dist[0]>3.5 or dist[1]>3.5 or dist[2]>3.5: continue

                                        dist = np.dot(dist,dist)
                                        if at_n[i] != 'H' and at_n[j] != 'H':
                                                #maximum bond length: 1.9A
                                                if 0.57 < dist < 12.25:
                                                        self._pbc_bonds[k].append([at_c[i]+l[0],at_c[j]+l[1],at_n[i],at_n[i]])
                                        else:
                                                #maximum bond length for hydrogen: 1.2A
                                                if 0.57 < dist < 5.15:
                                                        self._pbc_bonds[k].append([at_c[i]+l[0],at_c[j]+l[1],at_n[i],at_n[i]])

        #####################################################
        # EDIT FUNCTIONS
        #####################################################

        def mult(self,x,y,z):
                nat = self.get_nat()
                vec = self.get_vec()*self.get_celldm()
                mult = [x,y,z]
                for k in [0,1,2]:
                        for i in range(1,mult[k]):
                                for j in range(nat):
                                        self.append_atom_cp(self.get_atom(j))
                                        atom = self.get_atom(-1)
                                        self.set_atom(-1,atom[0],atom[1]+i*vec[k],'bohr')
                        nat = self.get_nat()
                self.set_vec(self.get_vec()*[[x],[y],[z]])

        #####################################################
        # VOLUME DATA FUNCTIONS
        #####################################################

        def set_vol(self,dim,vol):
                self.volume=np.array([[[0.]*dim[0]]*dim[1]]*dim[2],'f')
                i=0
                j=0
                line=vol[i].split()
                for x in range(dim[0]):
                        for y in range(dim[1]):
                                for z in range(dim[2]):
                                        self.volume[x][y][z]=float(line[j])
                                        j+=1
                                        if j==len(line) and i<(len(vol)-1):
                                                j=0
                                                i+=1
                                                line=vol[i].split()

        def vol_plane(self,height):
                if not (0<=height<self.nvol[2] or hasattr(self,'volume')):
                        return
                plane=self.volume[0:self.nvol[0],0:self.nvol[1],height]
                return plane

        #####################################################
        # SCRIPTING SUPPORT
        #####################################################

        def evalScript(self,script):
            script=script.split()
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
                    op=ops[script[0][0:3].lower()](script[1:])
                    #op=ops[script.pop(0)[0:3].lower()](script)
                except KeyError as e:
                    print 'Wrong Op:', e.message
                    return 'Wrong Op'
                except ValueError as e:
                    print 'Not an Angle!', e.message
                    return 'Not an Angle!'
                except NameError as e:
                    print 'Argument(s) missing, new command', e.message
                    return 'Argument(s) missing, new command', e.message
                except IndexError as e:
                    print 'Argument(s) missing, script over', e.message
                    return 'Argument(s) missing'
                except TypeError as e:
                    print e.message
                    return e.message
                else:
                    if op[0]=='define':
                        self._script_group[op[2]]=op[1]
                    elif op[0]=='repeat':
                        rep=op[1]
                    else:
                        for i in range(rep):
                            stack.append((op))
                        rep=1
                    del script[0:len(op)]
            # delete previous definitions
            self._script_group={}
            #if everything went well, execute operations
            for op in stack:
                op[0](*op[1:])
            return 'Success!'

        def _evalArgs(self,op,args):
            def evArgs(arglist):
                res = [op]
                for i,t in enumerate(args):
                    # if there's an argument, evaluate it
                    if i < len(arglist):
                        arg=arglist[i]
                    else:
                        #if vec is optional, command can be executed anyways
                        if t=='o':
                            return res
                        #else, there's an error
                        else:
                            raise IndexError('list index out of range')
                    #list of atoms
                    if t == 'l':
                        if arg in self._script_group:
                                res.append(self._script_group[arg])
                        elif type(eval(arg))==type(1):
                                res.append([int(arg)])
                        elif type(eval(arg))==type([]):
                                res.append(map(int,eval(arg)))
                        else:
                                raise TypeError('Not a list of atoms: '+str(arg))
                    #angle
                    elif t == 'a':
                        res.append(float(arg))
                    #index for loops
                    elif t == 'i':
                        res.append(int(arg))
                    #arbitrary names for defined groups
                    elif t == 's':
                        res.append(arg)
                    #valid vector for operations
                    elif t in 'vo':
                        arg=eval(arg)
                        if type(arg)==type(1) or (type(arg)==type([]) and len(arg)==1):
                                arg=self._atom_coord[arg-1]
                        elif type(arg)==type([]) and len(arg)==2:
                                arg=self._atom_coord[arg[0]]-self._atom_coord[arg[1]]
                        elif type(arg)==type(()) and len(arg)==4:
                                arg=self._set_coord(arg[0:3],arg[3])
                        elif type(arg)==type(()) and len(arg)==3:
                                arg=self._set_coord(arg)
                        else:
                                raise TypeError('Not a valid vector: '+str(arg))
                        res.append(arg)
                return res
            return evArgs

        def _rotate(self,atoms,angle,ax,shift=np.zeros(3)):
            angle=np.radians(angle)
            c=np.float(np.cos(angle))
            s=np.float(-np.sin(angle))
            ic=np.float(1.-c)
            mat=np.array([[ic*ax[0]*ax[0]+c,ic*ax[0]*ax[1]-s*ax[2],ic*ax[0]*ax[2]+s*ax[1]],
                          [ic*ax[0]*ax[1]+s*ax[2],ic*ax[1]*ax[1]+c,ic*ax[1]*ax[2]-s*ax[0]],
                          [ic*ax[0]*ax[2]-s*ax[1],ic*ax[1]*ax[2]+s*ax[0],ic*ax[2]*ax[2]+c]],'f')
            for i in atoms:
                self._atom_coord[i-1]=np.dot(self._atom_coord[i-1]-shift,mat)+shift

        def _shift(self,atoms,vector):
            for i in atoms:
                self._atom_coord[i-1]+=np.array(vector,'f')

        def _mirror(self,atoms,v1,v2,shift=np.zeros(3)):
            normal=np.cross(v1,v2)
            normal=normal/np.linalg.norm(normal)
            for i in atoms:
                pos=self._atom_coord[i-1]
                proj = np.dot(pos-shift,normal)*normal
                self._atom_coord[i-1]=pos-2*proj

