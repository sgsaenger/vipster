#!/usr/bin/env python
##########################################################################
# Molecule Class
# Atom Class
##########################################################################
version=0.3
versiontext='# class_molecule.py version {:.1f}'.format(version)
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
import class_molecule as cm

######################################################################
# MOLECULE CLASS
######################################################################
class molecule(cm.molecule):
            
    # create nanocrystal
    def crystal_make(self, mol,mx,my,mz,basis=[[0.0,0.0,0.0]],
                rota=0.0,rotv=[1.0,0.0,0.0],rotp=[0.0,0.0,0.0]):
        # copy the molecules
        v=self.Rvec()
        for ix in range(-mx,mx):
            for iy in range(-my,my):
                for iz in range(-mz,mz):
                    for ib in range(0,len(basis)):
                        x=float(ix)+float(basis[ib][0])
                        y=float(iy)+float(basis[ib][1])
                        z=float(iz)+float(basis[ib][2])
                        shift=[
                            x*v[0][0]+y*v[1][0]+z*v[2][0],
                            x*v[0][1]+y*v[1][1]+z*v[2][1],
                            x*v[0][2]+y*v[1][2]+z*v[2][2]]
                        self.append_submol(mol,shift,rota,rotv,rotp)
                        self.mol[self.Rnmol()-1].basis=ib
        return
     
    # cut crystal
    def crystal_cut(self,point,normal): # plane=a*x+b*y+c*z+f=0
        poplist=[]
        #f =(-ax0-by0-cz0)
        normal=calc.norm(normal)
        f=(-calc.vecprod(normal,point))
        # loop over molecules
        for i in range(0,self.Rnmol()):
            for j in range(0,self.mol[i].Rnatoms()): #len(self.mol[i].at)
                coo=self.mol[i].at[j].Rcoord()
                D=(calc.vecprod(normal,coo)+f)/calc.length(normal)
                # if above plane remove
                if D>0.0:
                    poplist.append(self.mol[i].Rid())
        # rm double or tripple counts
        poplist=sorted(list( set(poplist) ))
        # pop molecules
        for popi in range(0,len(poplist)):
            self.rm_submol(poplist[popi])            
        print ('... {:5d} molecules popped').format(len(poplist))
        return

