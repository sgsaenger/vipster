#!/usr/bin/env python
##########################################################################
# Molecule Class
# subclasses:
#   Atom Class
#   Bond Class
##########################################################################
version=1.7
versiontext='# mod_extend.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import copy
from operator import itemgetter, attrgetter
import mod_calc as calc # several functions

#----------------------------------------------------------------------
# classes
#----------------------------------------------------------------------


######################################################################
# MOLECULE CLASS
######################################################################
class molecule_extend():
    

    # extend molecule
    def extend(self):
        if not hasattr(self, "mol"):
            self.mol=[]
            # shift data
            self.mol.append(self.__class__())
            self.mol[0].at            = self.at
            self.mol[0].id            = self.id
            self.mol[0].natoms        = self.natoms
            self.mol[0].ntypes        = self.ntypes
            self.mol[0].typelist      = self.typelist
            self.mol[0].file          = self.file
            self.mol[0].filemolnumber = self.filemolnumber
            self.mol[0].vec           = self.vec
            # delete data
            self.file=""
            self.filenumber=0
            self.comment=0
            self.at=[]
            # print info
            self.extend_set()
            print >> sys.stderr, ('... molecule extended').format()
        else:
            print >> sys.stderr, ('... molecule already extended').format()
   
    # set extended molecule
    def extend_set(self):
        if hasattr(self, "mol"):
            # delete list
            for i in range(0,len(self.at)):
                self.at.pop()
            # make new list
            for i in range(0,len(self.at)):
                self.at.append(self.at[i])
            for i in range(0,len(self.mol)):
                for j in range(0,len(self.mol[i].at)):
                    self.at.append(self.mol[i].at[j])
            # calc natoms
            self.natoms=len(self.at)
            return
        else:
            return

    # append submolecule
    def append_submol(self,molecule,shiftv=[0.0,0.0,0.0],
                rotangle=0.0,rota=[1.0,0.0,0.0],rotp=[0.0,0.0,0.0]):
        # molecule information
        molid=len(self.mol)
        if hasattr(self,"mol"):
            # rotate and shift molecules
            molecule.rot(rotangle,rota[0],rota[1],rota[2],rotp[0],rotp[1],rotp[2])
            molecule.shift(shiftv[0],shiftv[1],shiftv[2])
            # create new submolecule and copy data
            self.mol.append(self.__class__())
            self.mol[molid].id=molid        # ID
            self.mol[molid].shiftvec=shiftv # shiftvec
            self.mol[molid].vec=copy.deepcopy(molecule.vec)
            # copy atoms
            for i in range (0,len(molecule.at)):
                at=molecule.at[i]
                self.mol[molid].at.append(
                    self.__class__.atom(self.mol[molid],i,at.name,at.number,
                                        at.coord[0],at.coord[1],at.coord[2],
                                        at.mult[0],at.mult[1],at.mult[2],
                                        at.charge)
                    )
            # set new submolecule
            self.mol[molid].set()
            # shift and rotate back     
            molecule.shift(-shiftv[0],-shiftv[1],-shiftv[2])
            molecule.rot(-rotangle,rota[0],rota[1],rota[2],rotp[0],rotp[1],rotp[2])
            # final set
            self.extend_set()
        else:
            print >> sys.stderr, ('... can not append submolecule').format()
        return

    # multiply molecule
    def mol_multiply(self,mx,my,mz):
        if not hasattr(self,"mol"): self.extend()
        # copy the molecules
        v=self.vec
        for ix in range(0,mx):
            for iy in range(0,my):
                for iz in range(0,mz):
                    x=float(ix)
                    y=float(iy)
                    z=float(iz)
                    shift=[x*v[0][0]+y*v[1][0]+z*v[2][0],
                           x*v[0][1]+y*v[1][1]+z*v[2][1],
                           x*v[0][2]+y*v[1][2]+z*v[2][2]]
                    if not(x==0 and y==0 and z==0): self.append_submol(self.mol[0],shift)
        return

    #remove submolecules
    def rm_submol(self,del_id):
        for i in range(0, len(self.mol)):
            if del_id==self.mol[i].id:
                self.mol.pop(i)
                break
        self.extend_set()
        return

    # rot all submolecules in list
    def mol_rot(self,angle,ax,ay,az,px,py,pz,lst=[]):
        if lst==[]:
            lst=xrange(1,len(self.mol))
        for imol in range(1,len(self.mol)):
            if self.mol[imol].id in lst:
                s=self.mol[imol].shiftvec
                self.mol[imol].rot(angle,ax,ay,az,px+s[0],py+s[1],pz+s[2])
        return

    # subgroup atoms
    def atoms_into_submol(self,idlist):
        if hasattr(self,"mol"):
            molid=len(self.mol)
            self.mol.append(self.__class__())
            # search for all atoms in idlist
            for idcnt in range(0,len(idlist)):
                found=False
                # search through all atoms
                for iat in range(0,len(self.mol[0].at)):
                    # if id equal get it
                    m=self.mol
                    if m[0].at[iat].id == idlist[idcnt]:
                        m[molid].append_atom_cp(m[0].at[iat])
                        m[molid].at[len(m[molid].at)-1].mol=m[molid]
                        # pop atom list
                        found=True
                        self.mol[0].at.pop(iat)
                        break
                # check if atom found
                if not found:
                    print >> sys.stderr, ('...atoms not found').format()
                    exit()
            # set shiftvec in molecule
            self.mol[molid].shiftvec=self.mol[molid].at[0].coord
        else:
            print >> sys.stderr, ('...molecule has to be extended').format()
            exit()
        return
