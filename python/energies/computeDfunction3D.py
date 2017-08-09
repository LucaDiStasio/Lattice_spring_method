# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def computeDfunction3D(t=None,N=None,lattice=None,dissipationparam=None,vel=None,dsthandle=None,dshhandle=None,dbehandle=None,structuralneighbours=None,shearneighbours=None,bendneighbours=None,*args,**kwargs):
    varargin = computeDfunction3D.varargin
    nargin = computeDfunction3D.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: July 2nd, 2014
#    Last update: August 4th, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    Dfuncp=zeros(N,1)
    # ---> structural springs
    
    for i in arange(1,N).reshape(-1):
        for j in arange(1,6).reshape(-1):
            neigh=structuralneighbours[i,j + 1]
            if neigh != - 1:
                Dfuncp[i,:]=Dfuncp[i,:] + multiply(dot(0.5,dsthandle[dissipationparam,t,dot(0.5,(lattice[neigh,7:9] + lattice[i,7:9])),lattice]),sum((multiply(sum((multiply((vel[i,:] - vel[neigh,:]),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))),2),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))) ** 2,2))
    
    # ---> shear springs
    
    for i in arange(1,N).reshape(-1):
        for j in arange(1,20).reshape(-1):
            neigh=shearneighbours[i,j + 1]
            if neigh != - 1:
                Dfuncp[i,:]=Dfuncp[i,:] + multiply(dot(0.5,dshhandle[dissipationparam,t,dot(0.5,(lattice[neigh,7:9] + lattice[i,7:9])),lattice]),sum((multiply(sum((multiply((vel[i,:] - vel[neigh,:]),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))),2),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))) ** 2,2))
    
    # ---> bend springs
    
    for i in arange(1,N).reshape(-1):
        for j in arange(1,6).reshape(-1):
            neigh=bendneighbours[i,j + 1]
            if neigh != - 1:
                Dfuncp[i,:]=Dfuncp[i,:] + multiply(dot(0.5,dbehandle[dissipationparam,t,dot(0.5,(lattice[neigh,7:9] + lattice[i,7:9])),lattice]),sum((multiply(sum((multiply((vel[i,:] - vel[neigh,:]),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))),2),(((lattice[neigh,7:9] - lattice[i,7:9])) / (sqrt(sum((lattice[neigh,7:9] - lattice[i,7:9]) ** 2,2)))))) ** 2,2))
    
    Dfunc=sum(Dfuncp)
    return Dfuncp,Dfunc