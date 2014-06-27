#!/usr/bin/env python
##########################################################################
# writepwscf
# readfpwscfin
# readfpwscfout
##########################################################################
version=1.0
versiontext='# mod_pwscf.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import string
import copy
import mod_calc as calc
from numpy import matrix
from numpy import linalg

#----------------------------------------------------------------------
# module mod_xyz
#----------------------------------------------------------------------
class molecule_rw:
    #############################################################
    # return functions
    #############################################################                      
    def __pwscfinit__(self):
        self.celldm=float(1.0)
    
    def Rcelldm(self):
        return self.celldm

    #############################################################
    # modify functions
    #############################################################    
    # set celldm
    def set_celldm(self,celldm):
        self.celldm=celldm
        return

    # write pwscf output format
    def writepwscf(self,filename="",status='w'):
        # calculate relative coordinates
        if not hasattr(self.at[0],"coord_rel"):
            vecM = matrix( [ [self.Rvec()[0][0],self.Rvec()[1][0],self.Rvec()[2][0]],
                             [self.Rvec()[0][1],self.Rvec()[1][1],self.Rvec()[2][1]],
                             [self.Rvec()[0][2],self.Rvec()[1][2],self.Rvec()[2][2]] ])
            for cntat in range(self.Rnatoms()):
                tmp = matrix([ [self.at[cntat].Rcoord()[0]],
                               [self.at[cntat].Rcoord()[1]],
                               [self.at[cntat].Rcoord()[2]] ])
                res=linalg.solve(vecM, tmp)
                self.at[cntat].coord_rel=[float(res[0]),float(res[1]),float(res[2])]
                print self.at[cntat].coord_rel
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
        mol=self
        # print options if present
        if hasattr(mol, 'control'):
            print >>f, "&control"
            for i in range(len(mol.control)):
                print >>f, ('{:s}={:s}').format(mol.control[i][0],mol.control[i][1])
            print >>f, "/"
        if hasattr(mol, 'system'):
            print >>f, "&system"
            for i in range(len(mol.system)):
                if mol.system[i][0]=="celldm(1)": 
                    print >>f, ('{:s}={:f}').format(mol.system[i][0],mol.Rcelldm()/calc.b2A)
                else:
                    print >>f, ('{:s}={:s}').format(mol.system[i][0],mol.system[i][1])
            print >>f, "/"
        if hasattr(mol, 'electrons'):
            print >>f, "&electrons"
            for i in range(len(mol.electrons)):
                print >>f, ('{:s}={:s}').format(mol.electrons[i][0],mol.electrons[i][1])
            print >>f, "/"
        if hasattr(mol, 'ions'):
            print >>f, "&ions"
            for i in range(len(mol.ions)):
                print >>f, ('{:s}={:s}').format(mol.ions[i][0],mol.ions[i][1])
            print >>f, "/"
        # print species if present
        # TODO otherwise calculate them at the beginning
        if hasattr(mol, 'species'):
            print >>f
            print >>f, "ATOMIC_SPECIES crystal"
            for i in range(len(mol.species)):
                print >>f, ('{:s}').format(mol.species[i]),
            print >>f
        # print relative coordinates
        print >>f
        print >>f, "ATOMIC_POSITIONS"
        for cntat in range(mol.Rnatoms()):
            print >>f ,(
                '{:4s} {:15.10f} {:15.10f} {:15.10f}'.format(
                    mol.at[cntat].Rname(),
                    mol.at[cntat].coord_rel[0], 
                    mol.at[cntat].coord_rel[1], 
                    mol.at[cntat].coord_rel[2]
                    )
                )
        print >>f
        print >>f, "CELL_PARAMETERS"
        for cntvec in range(0,3):
            print >>f , (
                '{:15.10f} {:15.10f} {:15.10f}'.format(
                    mol.Rvec()[cntvec][0]/mol.Rcelldm(), 
                    mol.Rvec()[cntvec][1]/mol.Rcelldm(), 
                    mol.Rvec()[cntvec][2]/mol.Rcelldm()
                    )
                )
        if hasattr(mol,"kpoints"):
            print >>f
            print >>f, "K_POINTS"
            print >>f , ('{:s}'.format(mol.kpoints[0]))
        f.close()
        return
    
    # read molecules in pwscf input file format
    def readpwscfin(self,filename):
        # set molecule
        molecules=[]
        natoms=0
        vec=[]
        # read file
        file=open(filename, 'r')
        opt=""
        cntline=0
        cntat=0
        cntvec=0
        cnttypes=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            #
            # create new molecule
            #
            if cntline==1:
                mol=self.__class__()
                mol.at=[]
                # additional optional fields
                mol.control=[]
                mol.system=[]
                mol.electrons=[]
                mol.ions=[]
                mol.kpoints=[]
                mol.species=[]
            #
            # Do READ IN
            #
            # read vectors
            if  opt=="readvec":
                vec.append([float(linesplit[0])*mol.Rcelldm(),
                            float(linesplit[1])*mol.Rcelldm(),
                            float(linesplit[2])*mol.Rcelldm()])
                cntvec+=1
                if cntvec==3: opt=""
            # read coordinates
            if  opt=="readspecies":
                mol.species.append(line)
                cnttypes+=1
                if cnttypes==ntypes: opt=""
            # read unitvector
            elif opt=="readcoord":
                mol.at.append(
                    self.__class__.atom(
                        mol,
                        cntat,
                        linesplit[0],0,
                        float(linesplit[1]),
                        float(linesplit[2]),
                        float(linesplit[3])
                        )
                    )
                cntat+=1
                if int(cntat)==int(natoms): opt=""
            # read input options
            elif opt=="system":
                if linesplit[0]=="/": 
                    opt=""
                    for i in range(len(mol.system)):
                        s=mol.system[i]
                        if   s[0]=="nat":       natoms=int(s[1])
                        elif s[0]=="ntyp":      ntypes=int(s[1])
                        elif s[0]=="celldm(1)": mol.set_celldm(float(s[1])*calc.b2A)
                else:
                    for i in range(len(linesplit)):
                        linesplit[i]=linesplit[i].replace(",","")
                        mol.system.append(linesplit[i].split("="))
            elif opt=="control":
                if linesplit[0]=="/": opt=""
                else:
                    for i in range(len(linesplit)):
                        linesplit[i]=linesplit[i].replace(",","")
                        mol.control.append(linesplit[i].split("="))
            elif opt=="electrons":
                if linesplit[0]=="/": opt=""
                else:
                    for i in range(len(linesplit)):
                        linesplit[i]=linesplit[i].replace(",","")
                        mol.electrons.append(linesplit[i].split("="))
            elif opt=="ions":
                if linesplit[0]=="/": opt=""
                else:
                    for i in range(len(linesplit)):
                        linesplit[i]=linesplit[i].replace(",","")
                        mol.ions.append(linesplit[i].split("="))
            elif opt=="readkpoints":
                mol.kpoints.append(line)
                opt=""
            #
            # read main options
            #
            if len(linesplit)>0:
                option=linesplit[0]
                if  option=="&control":
                    opt="control"
                elif option=="&system":
                    opt="system"
                elif option=="&electrons":
                    opt="electrons"
                elif option=="&ions":
                    opt="ions"
                elif option=="ATOMIC_SPECIES":
                    opt="readspecies"
                elif option=="ATOMIC_POSITIONS":
                    opt="readcoord"
                elif option=="CELL_PARAMETERS":
                    opt="readvec"
                elif option=="K_POINTS":
                    opt="readkpoints"                    
        # set periodicity
        mol.set_periodicity(vec[0],vec[1],vec[2])
        # set real coordinates
        mol.rel2real()
        # set molecule and append
        mol.set(filename,1,"")
        molecules.append(copy.copy(mol))
        # close file
        file.close()
        # return molecules
        return molecules

    # read molecules in pwscf output file format
    def readpwscfout(self,filename):
        # set molecule
        molecules=[]
        natoms=0
        vec=[]
        # read file
        file=open(filename, 'r')
        opt=""
        cntline=0
        oldline=0
        cntat=0
        cntvec=0
        cnttypes=0
        for line in file:
            cntline+=1
            linesplit=line.split()
            #
            # create new molecule
            #
            if  cntline==1:
                mol=self.__class__()
                mol.at=[]
                # additional optional fields
                #mol.control=[]
                #mol.system=[]
                #mol.electrons=[]
                #mol.ions=[]
                #mol.kpoints=[]
                mol.species=[]
                    
            if  opt=="readcoord" and cntat==0:
                if not len(mol.at)==0:
                    mol=self.__class__()
                    mol.at=[]
                    # additional optional fields
                    #mol.control=[]
                    #mol.system=[]
                    #mol.electrons=[]
                    #mol.ions=[]
                    #mol.kpoints=[]
                    #mol.species=[]
                # attribute global attributes to molecule
                mol.set_celldm(celldm*calc.b2A)
            #
            # Do READ IN
            #
            # read vectors
            if  opt=="readvec":
                vec.append([float(linesplit[3]),
                            float(linesplit[4]),
                            float(linesplit[5])])
                cntvec+=1
                if cntvec==3: opt=""
            # read coordinates
            if  opt=="readspecies":
                cnttypes+=1
                if cnttypes==ntypes: opt=""
            # read unitvector
            elif opt=="readcoord":
                mol.at.append(
                    self.__class__.atom(
                        mol,
                        cntat,
                        linesplit[0],0,
                        float(linesplit[1]),
                        float(linesplit[2]),
                        float(linesplit[3])
                        )
                    )
                cntat+=1
                if int(cntat)==int(natoms): 
                    opt=""
                    cntat=0
                    # append mol to molecules
                    # set periodicity
                    mol.set_periodicity(calc.scal_vecmult(mol.Rcelldm(),vec[0]) ,
                                        calc.scal_vecmult(mol.Rcelldm(),vec[1]) ,
                                        calc.scal_vecmult(mol.Rcelldm(),vec[2]) )
                    # set real coordinates
                    mol.rel2real()
                    # set molecule and append
                    mol.set(filename,1,"")
                    molecules.append(copy.copy(mol))
            # read input options
            elif opt=="":
                if   len(linesplit)>1 and linesplit[0:3]==["number","of","atoms/cell"]:
                    natoms=int(linesplit[4])
                elif len(linesplit)>3 and linesplit[0:4]==["number","of","atomic","types"]:
                    ntypes=int(linesplit[5])                        
                elif len(linesplit)>1 and linesplit[0]=="celldm(1)=":
                    mol.set_celldm(float(linesplit[1])*calc.b2A)
            #
            # read main options
            #
            if len(linesplit)>0:
                option=linesplit[0]
                if   option=="atomic" and linesplit[1]=="species":
                    opt="readspecies"
                elif option=="ATOMIC_POSITIONS":
                    opt="readcoord"
                elif option=="crystal" and linesplit[1]=="axes:":
                    opt="readvec"
        # close file
        file.close()
        # return molecules
        return molecules

    def rel2real(self):
        for i in range(len(self.at)):
            coo=[float(0.0),float(0.0),float(0.0)]
            # calculate coordinates
            for dim in range(len(self.at[i].coord)):
                coo=calc.vecadd(coo,calc.scal_vecmult(self.at[i].coord[dim],self.vec[dim]))
            # set relative and real coordinates
            self.at[i].coord_rel=self.at[i].coord
            self.at[i].spos(coo[0],coo[1],coo[2])
