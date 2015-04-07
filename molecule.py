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
                self._bond_cutoff=[]
                self._script_group=dict()
                self._celldm = 1.0
                self._vec=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],'f')
                self._vecinv=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]],'f')
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
        def set_celldm(self,cdm,scale=False):
                if scale:
                    ratio=cdm/self._celldm
                    for i in range(len(self._atom_coord)):
                        self._atom_coord[i] = self._atom_coord[i]*ratio
                self._celldm = float(cdm)

        # set vectors in bohr
        def set_vec(self,vec,scale=False):
                vec = np.array(vec,'f')
                inv = self._vecinv
                if scale:
                    for i in range(len(self._atom_coord)):
                        self._atom_coord[i] = np.dot(np.dot(self._atom_coord[i],inv),vec)
                self._vec = vec
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

        def get_bonds(self):
                if not hasattr(self,'_bonds'): self.set_bonds()
                return self._bonds

        def set_bonds(self):
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
                except KeyError as e:
                    return 'Wrong Op: '+script[0]
                except NameError as e:
                    print e.message
                    return e.message
                except TypeError as e:
                    print e.message
                    return e.message
                except IndexError as e:
                    print e.message
                    return e.message
                except ValueError as e:
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
                            raise IndexError('Argument index out of range')
                    #list of atoms
                    if t == 'l':
                        if arg in self._script_group:
                            res.append(self._script_group[arg])
                        else:
                            arg = arg.strip('[]').split(',')
                            l = []
                            for j in arg:
                                if '-' in j:
                                    low,high=j.split('-')
                                    l.extend(range(int(low)-1,int(high)))
                                else:
                                    l.append(int(j)-1)
                            if not np.all(np.less(l,len(self._atom_coord))):
                                raise IndexError('Atom index out of range')
                            res.append(l)
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
                        if '(' in arg:
                            arg=eval(arg)
                            if len(arg)==4:
                                arg=self._set_coord(arg[0:3],arg[3])
                            else:
                                arg=self._set_coord(arg[0:3])
                        elif '-' in arg:
                            arg=arg.strip('[]').split('-')
                            arg=self._atom_coord[int(arg[0])-1]-self._atom_coord[int(arg[1])-1]
                        elif type(eval(arg))==type(1) or (type(arg)==type([]) and len(arg)==1):
                            arg=self._atom_coord[int(arg)-1]
                        else:
                            raise ValueError('Not a valid vector: '+str(arg))
                        res.append(arg)
                return res
            return evArgs

        def _rotate(self,atoms,angle,ax,shift=np.zeros(3)):
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
            for i in atoms:
                self._atom_coord[i]+=np.array(vector,'f')

        def _mirror(self,atoms,v1,v2,shift=np.zeros(3)):
            normal=np.cross(v1,v2)
            normal=normal/np.linalg.norm(normal)
            for i in atoms:
                pos=self._atom_coord[i]
                proj = np.dot(pos-shift,normal)*normal
                self._atom_coord[i]=pos-2*proj

