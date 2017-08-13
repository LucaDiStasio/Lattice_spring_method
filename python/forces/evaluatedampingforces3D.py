# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def evaluatedampingforces3D(t=None,N=None,lattice=None,dissipationparam=None,vel=None,dsthandle=None,dshhandle=None,dbehandle=None,structuralneighbours=None,shearneighbours=None,bendneighbours=None,*args,**kwargs):
    varargin = evaluatedampingforces3D.varargin
    nargin = evaluatedampingforces3D.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: June 26th, 2014
#    Last update: August 4th, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    Fdamp=zeros(N,3)
    ##
    
    # ---> structural springs
    
    for i in arange(1,N).reshape(-1):
        for j in arange(1,6).reshape(-1):
            neigh=structuralneighbours[i,j + 1]
            if neigh != - 1:
                Fdamp[i,:]=Fdamp[i,:] + multiply(multiply(dsthandle[dissipationparam,t,dot(0.5,(lattice[neigh,7:9] + lattice[i,7:9])),lattice],sum((multiply((vel[i,:] - vel[neigh,:]),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))),2)),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))
    
    # ---> shear springs
    
    for i in arange(1,N).reshape(-1):
        for j in arange(1,20).reshape(-1):
            neigh=shearneighbours[i,j + 1]
            if neigh != - 1:
                Fdamp[i,:]=Fdamp[i,:] + multiply(multiply(dshhandle[dissipationparam,t,dot(0.5,(lattice[neigh,7:9] + lattice[i,7:9])),lattice],sum((multiply((vel[i,:] - vel[neigh,:]),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))),2)),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))
    
    # ---> bend springs
    
    for i in arange(1,N).reshape(-1):
        for j in arange(1,6).reshape(-1):
            neigh=bendneighbours[i,j + 1]
            if neigh != - 1:
                Fdamp[i,:]=Fdamp[i,:] + multiply(multiply(dbehandle[dissipationparam,t,dot(0.5,(lattice[neigh,7:9] + lattice[i,7:9])),lattice],sum((multiply((vel[i,:] - vel[neigh,:]),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))),2)),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))
    
    return Fdamp