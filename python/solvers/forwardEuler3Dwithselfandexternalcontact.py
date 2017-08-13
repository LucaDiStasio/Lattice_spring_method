# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def forwardEuler3Dwithselfandexternalcontact(elasticparam=None,dissipationparam=None,tn=None,dt=None,N=None,Nx=None,Ny=None,lattice=None,vel=None,m=None,ksthandle=None,kshhandle=None,kbehandle=None,dsthandle=None,dshhandle=None,dbehandle=None,fexthandle=None,boundarysets=None,boundarytype=None,posfuncs=None,velfuncs=None,indicesbulk=None,indicesinternalbulk=None,indicesF1=None,indicesF2=None,indicesF3=None,indicesF4=None,indicesF5=None,indicesF6=None,indicesinternalF1=None,indicesinternalF2=None,indicesinternalF3=None,indicesinternalF4=None,indicesinternalF5=None,indicesinternalF6=None,indicesE1=None,indicesE2=None,indicesE3=None,indicesE4=None,indicesE5=None,indicesE6=None,indicesE7=None,indicesE8=None,indicesE9=None,indicesE10=None,indicesE11=None,indicesE12=None,indicesinternalE1=None,indicesinternalE2=None,indicesinternalE3=None,indicesinternalE4=None,indicesinternalE5=None,indicesinternalE6=None,indicesinternalE7=None,indicesinternalE8=None,indicesinternalE9=None,indicesinternalE10=None,indicesinternalE11=None,indicesinternalE12=None,indicesC1=None,indicesC2=None,indicesC3=None,indicesC4=None,indicesC5=None,indicesC6=None,indicesC7=None,indicesC8=None,indicesinternalC1=None,indicesinternalC2=None,indicesinternalC3=None,indicesinternalC4=None,indicesinternalC5=None,indicesinternalC6=None,indicesinternalC7=None,indicesinternalC8=None,*args,**kwargs):
    varargin = forwardEuler3Dwithselfandexternalcontact.varargin
    nargin = forwardEuler3Dwithselfandexternalcontact.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: August 4th, 2014
#    Last update: August 4th, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    Fel=evaluateelasticforces3D(tn,N,Nx,Ny,lattice,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8)
    Fdamp=evaluatedampingforces3D(tn,N,Nx,Ny,lattice,dissipationparam,vel,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8)
    Fext=evaluateexternalforces3D(tn,lattice,vel,fexthandle)
    Fselfcontact=evaluateselfcontact3D(N,lattice,boundary,Rcsq,ksc)
    Fextcontact=evaluateexternalcontact3D(t,N,lattice,boundary,boundhandle,gradienthandle,kec)
    lattice[:,7:9]=lattice[:,7:9] + multiply(dt,vel)
    vel=vel + multiply(dt,(Fel + Fdamp + Fext + Fselfcontact + Fextcontact)) / m
    lattice,vel=set_boundaryconditions(tn + dt,N,lattice,vel,boundarysets,boundarytype,posfuncs,velfuncs,nargout=2)
    return lattice,vel