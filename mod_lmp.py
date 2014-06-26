#!/usr/bin/env python
##########################################################################
# writexyz
# readfxyz
##########################################################################
version=0.9
versiontext='# mod_lmp.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
import copy
import mod_calc as calc

#----------------------------------------------------------------------
# module mod_xyz
#----------------------------------------------------------------------
class molecule_rw:
    # write lammps file
    def writelmp(self,filename="",status='w'):
        ndim=calc.ndim
        # define vec and offset
        v=self.vec
        o=self.offset
        # check if vector is correct for lammps
        if not (v[1][0]==0.0 and v[2][0]==0.0 and v[2][1]==0.0): 
            print "vectors are not defined correctly for a lammps file"
            quit()
        # open file if present
        if filename == "":
            f=sys.stdout
        else:
            f=open(filename, status)
        mol=self
        print >>f
        print >>f
        # numbers
        print >>f, ("{:d} atoms".format(mol.natoms))
        print >>f
        # types
        print >>f, ("{:d} atom types".format(mol.ntypes))
        print >>f
        # vectors
        if not (v[0][1]==0.0 and v[0][2]==0.0 and v[1][2]==0.0):
            print >>f, ("{:f} {:f} {:f} xy xz yz".format(v[0][1],v[0][2],v[1][2]))
            print >>f
        # offset
        print >>f, ("{:20.10f} {:20.10f} xlo xhi".format(o[0],v[0][0]+o[0]))
        print >>f, ("{:20.10f} {:20.10f} ylo yhi".format(o[1],v[1][1]+o[1]))
        print >>f, ("{:20.10f} {:20.10f} zlo zhi".format(o[2],v[2][2]+o[2]))
        print >>f
        # atomtypes # CHANGE
        print >> f, "Masses"
        for cnttype in range(mol.ntypes):
            print >>f ,(
                '{:6d} {:f} #{:s}'.format(
                    cnttype+1,
                    self.type_number2weight(mol.typelist[cnttype]),
                    self.type_weight2name(self.type_number2weight(mol.typelist[cnttype])),
                    )
                )       
        print >>f
        # atoms
        print >> f, "Atoms"
        for cntat in range(mol.natoms):
            print >>f ,(
                '{:6d} {:4d} {:15.10f} {:15.10f} {:15.10f}'.format(
                    cntat+1, 
                    mol.at[cntat].tid+1,
                    mol.at[cntat].coord[0], 
                    mol.at[cntat].coord[1], 
                    mol.at[cntat].coord[2]
                    )
                )
        print >>f
        f.close()
        return
    
    # read molecules in lammps file
    def readlmp(self,filename):
        ndim=calc.ndim
        # set molecule
        molecules=[]
        tilt=[0.0 for j in range(ndim)]
        vec=[[0.0 for j in range(ndim)] for i in range(ndim)]
        local_types=[]
        # read file
        file=open(filename, 'r')
        cntline=0
        natoms=0
        opt=""
        for line in file:
            cntline+=1
            linesplit=line.split()
            # create new molecule
            if (cntline)==1:
                mol=self.__class__()
                mol.at=[]
                cntat=0
                cnttype=0
            if len(linesplit)>0 and linesplit[0][0]=="#":
                continue
            # read atoms
            if opt=="atoms":
                # read number and name
                tid    = int(linesplit[1])
                name   = self.local_number2name(tid,local_types)
                number = self.type_name2number(name)
                # append atoms
                mol.at.append(
                    #cm.molecule.atom(
                    self.__class__.atom(
                        mol,
                        cntat,
                        name,
                        number,
                        float(linesplit[2]),
                        float(linesplit[3]),
                        float(linesplit[4]),
                        tid=tid-1
                        )
                    )
                cntat+=1
                if cntat==natoms: opt=""
            # read types
            if opt=="masses":
                # append types
                local_types.append([int(linesplit[0]),float(linesplit[1])])
                # until cnttype==ntypes
                cnttype+=1
                if cnttype==ntypes: opt=""
            #
            # read main options
            #
            # opt for blocks
            if len(linesplit)>0:
                option=linesplit[0]
                if option=="Masses":
                    opt="masses"
                elif option=="Atoms":
                    opt="atoms"
                elif option=="Bonds":
                    opt="bonds"
            # NUMBERS
            if len(linesplit)>=2:                
                option=linesplit[1]
                if   option=="atoms":
                    natoms=int(linesplit[0])
                elif option=="bonds":
                    nbonds=int(linesplit[0])
                elif option=="angles":
                    nangles=int(linesplit[0])
                elif option=="dihredrals":
                    ndihredrals=int(linesplit[0])
                elif option=="impropers":
                    nimpropers=int(linesplit[0])
            # TYPES
            if len(linesplit)>=3:                
                option=linesplit[1:3]
                if   option==["atom","types"]:
                    ntypes=int(linesplit[0])
                elif option==["bond","types"]:
                    tbonds=int(linesplit[0])
                elif option==["angle","types"]:
                    tangle=int(linesplit[0])
                elif option=="dihredrals":
                    tdihredrals=int(linesplit[0])
                elif option=="impropers":
                    timpropers=int(linesplit[0])
            # Vectors
            if len(linesplit)>=4:
                option=linesplit[2:4]
                if   option==["xlo","xhi"]:
                    x=[float(linesplit[0]),float(linesplit[1])]
                elif option==["ylo","yhi"]:
                    y=[float(linesplit[0]),float(linesplit[1])]
                elif option==["zlo","zhi"]:
                    z=[float(linesplit[0]),float(linesplit[1])]
            if len(linesplit)==6:
                if linesplit[3:6]==["xy","xz","yz"]:
                    tilt=[float(linesplit[0]),float(linesplit[1]),float(linesplit[2])]
        # calculate vectors a,b,c out of the box and tilt
        vec[0][:]=[ x[1]-x[0], tilt[0],      tilt[1]      ]
        vec[1][:]=[ 0.0,       y[1]-y[0],    tilt[2]      ]
        vec[2][:]=[ 0.0,       0.0,          z[1]-z[0]    ]
        # finish file and append it
        # set periodicity
        mol.set_periodicity(vec[0],vec[1],vec[2],[x[0],y[0],z[0]])
        # set molecule and append
        mol.set(filename,1,"")
        molecules.append(copy.copy(mol))
        # close file
        file.close()
        # return molecules
        return molecules
    
    def local_number2name(self,number,local_types):
        # local_types[type][0] id
        # local_types[type][2] ~ weight
        t=local_types
        name=""
        weight=0.0
        for cnt in range(len(t)):
            if number==t[cnt][0]:
                weight=t[cnt][1]
        # get name to return
        return self.type_weight2name(weight)

