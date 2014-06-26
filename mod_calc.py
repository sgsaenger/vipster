#!/usr/bin/env python
##########################################################################
# mod_calc module
##########################################################################
version=1.10
versiontext='# mod_calc.py version {:.1f}'.format(version)
#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math
from operator import itemgetter, attrgetter

#----------------------------------------------------------------------
# functions
#----------------------------------------------------------------------
b2A=0.52917725
ndim=3

# calculate distance
def a_dist(at1, at2, per=[0,0,0],dim=3):
    add=[0.0,0.0,0.0]
    vp=at1.mol.vec
    perf=[float(per[0]),float(per[1]),float(per[2])]
    add[0] = vp[0][0]*perf[0] + vp[1][0]*perf[1] + vp[2][0]*perf[2]
    add[1] = vp[0][1]*perf[0] + vp[1][1]*perf[1] + vp[2][1]*perf[2]
    add[2] = vp[0][2]*perf[0] + vp[1][2]*perf[1] + vp[2][2]*perf[2]
    d=[0.0 for x in range(0,dim)]
    for i in range(0,dim): d[i]=(at2.coord[i]-at1.coord[i]+add[i])
    dxyz = 0.0
    for i in range(0,dim): dxyz = dxyz + d[i]*d[i]
    return math.sqrt(dxyz)

def a_vec(at1,at2,per=[0,0,0]):
    add=[0.0,0.0,0.0]
    vp=at1.mol.vec
    perf=[float(per[0]),float(per[1]),float(per[2])]
    add[0] = vp[0][0]*perf[0] + vp[1][0]*perf[1] + vp[2][0]*perf[2]
    add[1] = vp[0][1]*perf[0] + vp[1][1]*perf[1] + vp[2][1]*perf[2]
    add[2] = vp[0][2]*perf[0] + vp[1][2]*perf[1] + vp[2][2]*perf[2]
    vec=[at2.coord[0]-at1.coord[0]+add[0],
         at2.coord[1]-at1.coord[1]+add[1],
         at2.coord[2]-at1.coord[2]+add[2]]
    return vec

def length(vec,dim=3):
    vecsum=0.0
    for cntdim in range(0,dim):
        vecsum= vecsum + (vec[cntdim]*vec[cntdim])
    vecsum=math.sqrt(vecsum)
    return vecsum

def vecsub(vec1,vec2,dim=3):
    vec=[0.0 for x in range(0,dim)]
    for cntdim in range(0,dim):
        vec[cntdim]=vec2[cntdim]-vec1[cntdim]
    return vec

def vecadd(vec1,vec2,dim=3):
    vec=[0.0 for x in range(0,dim)]
    for cntdim in range(0,dim):
        vec[cntdim]=vec2[cntdim]+vec1[cntdim]
    return vec           

def vecprod(vec1,vec2,dim=3):
    vecsum=0.0
    for cntdim in range(0,dim):
        vecsum+=vec2[cntdim]*vec1[cntdim]
    return vecsum

def vec_vecmult(vec1,vec2,dim=3):
    vec=[0.0 for x in range(0,dim)]
    for cntdim in range(0,dim):
        vec[cntdim]=vec1[cntdim]*vec2[cntdim]
    return vec           

def scal_vecmult(scal1,vec2,dim=3):
    vec=[0.0 for x in range(0,dim)]
    for cntdim in range(0,dim):
        vec[cntdim]=scal1*vec2[cntdim]
    return vec           

def norm(vec,dim=3):
    norm=[0.0 for x in xrange(0,dim)]
    l=length(vec,dim)
    for cntdim in range(0,dim):
        norm[cntdim]=vec[cntdim]/l
    return norm

def calc_crossproduct3d(vec1,vec2):
    tvec=[0.0,0.0,0.0]
    tvec[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1]
    tvec[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2]
    tvec[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0]
    return tvec

def mat2d_inv(mat):
    det= mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0]
    inv=[[0.0 for cnt in range(0,2)] for cnt in range(0,2)]
    inv[0][0] =  mat[1][1]/det
    inv[1][0] = -mat[0][1]/det
    inv[0][1] = -mat[1][0]/det
    inv[1][1] =  mat[0][0]/det
    return inv
           
def mat2d_mult(mat1,mat2):
    res = [[0.0 for i in range(0,2)] for i in range(0,2)]
    for dim1 in range(0,2):
        for dim2 in range(0,2):
            for cnt in range(0,2):
                res[dim1][dim2] = res[dim1][dim2] + (
                    mat1[cnt][dim2]*mat2[dim1][cnt]
                               )
    return res

def mat2d_transpose(mat):
    res = [[0.0 for i in range(0,2)] for i in range(0,2)]
    res[0][0] = mat[0][0]
    res[1][0] = mat[0][1]
    res[0][1] = mat[1][0]
    res[1][1] = mat[1][1]
    return res

def transform_vec(vec,unitvec,dim=3):
    # define unitvectors as
    #  ex=unitvec[0]
    #  ey=unitvec[1]
    #  ez=unitvec[2]   
    u=unitvec
    res = [0.0 for x in range(0,dim)]
    # do transformation
    for cnt1 in range(0,dim):
        for cnt2 in range(0,dim):
            res[cnt1] = res[cnt1] + u[cnt1][cnt2]*vec[cnt2]
    # return result
    return res

def data_deriv(data,delta=0.0):
    # calculate derivative of x:[0] and y:[1] of data
    # derivative is delta_y/delta_x at pos x_1+1/2*delta_x
    deriv=[]
    for i in range(0,len(data)-1):
        delta_x=(data[i+1][0]-data[i][0])
        if delta_x>delta:
            deriv.append([
                    (data[i][0]+0.5*delta_x)  ,
                    (data[i+1][1]-data[i][1])/delta_x
                    ])
    return deriv

def data_2nd_deriv(data):
    # calculate second derivative of x:[0] and y:[1] of data
    return data_deriv( data_deriv(data) )
