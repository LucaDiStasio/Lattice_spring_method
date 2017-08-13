# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def RK42D(elasticparam=None,dissipationparam=None,tn=None,dt=None,N=None,Nx=None,Ny=None,lattice=None,vel=None,m=None,ksthandle=None,kshhandle=None,kbehandle=None,dsthandle=None,dshhandle=None,dbehandle=None,fexthandle=None,boundarysets=None,boundarytype=None,posfuncs=None,velfuncs=None,indicesbulk=None,indicesinternalbulk=None,indicesF1=None,indicesF2=None,indicesF3=None,indicesF4=None,indicesF5=None,indicesF6=None,indicesinternalF1=None,indicesinternalF2=None,indicesinternalF3=None,indicesinternalF4=None,indicesinternalF5=None,indicesinternalF6=None,indicesE1=None,indicesE2=None,indicesE3=None,indicesE4=None,indicesE5=None,indicesE6=None,indicesE7=None,indicesE8=None,indicesE9=None,indicesE10=None,indicesE11=None,indicesE12=None,indicesinternalE1=None,indicesinternalE2=None,indicesinternalE3=None,indicesinternalE4=None,indicesinternalE5=None,indicesinternalE6=None,indicesinternalE7=None,indicesinternalE8=None,indicesinternalE9=None,indicesinternalE10=None,indicesinternalE11=None,indicesinternalE12=None,indicesC1=None,indicesC2=None,indicesC3=None,indicesC4=None,indicesC5=None,indicesC6=None,indicesC7=None,indicesC8=None,indicesinternalC1=None,indicesinternalC2=None,indicesinternalC3=None,indicesinternalC4=None,indicesinternalC5=None,indicesinternalC6=None,indicesinternalC7=None,indicesinternalC8=None,*args,**kwargs):
    varargin = RK42D.varargin
    nargin = RK42D.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: July 14th, 2014
#    Last update: August 4th, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    Fn=evaluateelasticforces2D(tn,N,Nx,Ny,lattice,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluatedampingforces2D(tn,N,Nx,Ny,lattice,dissipationparam,vel,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluateexternalforces2D(tn,lattice,vel,fexthandle)
    x1=lattice[:,5:6] + multiply(dot(0.5,dt),vel)
    vel1=vel + multiply(dot(0.5,dt),F) / m
    lattice1,vel1=set_boundaryconditions(tn + dt,N,cat(lattice[:,1:4],x1),vel1,boundarysets,boundarytype,posfuncs,velfuncs,nargout=2)
    F1=evaluateelasticforces2D(tn + dot(0.5,dt),N,Nx,Ny,lattice1,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluatedampingforces2D(tn + dot(0.5,dt),N,Nx,Ny,lattice1,dissipationparam,vel1,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluateexternalforces2D(tn + dot(0.5,dt),lattice1,vel1,fexthandle)
    x2=lattice[:,5:6] + multiply(dot(0.5,dt),vel1)
    vel2=vel + multiply(dot(0.5,dt),F1) / m
    lattice2,vel2=set_boundaryconditions(tn + dt,N,cat(lattice[:,1:4],x2),vel2,boundarysets,boundarytype,posfuncs,velfuncs,nargout=2)
    F2=evaluateelasticforces2D(tn + dot(0.5,dt),N,Nx,Ny,lattice2,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluatedampingforces2D(tn + dot(0.5,dt),N,Nx,Ny,lattice2,dissipationparam,vel2,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluateexternalforces2D(tn + dot(0.5,dt),lattice2,vel2,fexthandle)
    x3=lattice[:,5:6] + multiply(dt,vel2)
    vel3=vel + multiply(dt,F2) / m
    lattice3,vel3=set_boundaryconditions(tn + dt,N,cat(lattice[:,1:4],x3),vel3,boundarysets,boundarytype,posfuncs,velfuncs,nargout=2)
    F3=evaluateelasticforces2D(tn + dt,N,Nx,Ny,lattice3,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluatedampingforces2D(tn + dt,N,Nx,Ny,lattice3,dissipationparam,vel3,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8) + evaluateexternalforces2D(tn + dt,lattice3,vel3,fexthandle)
    lattice[:,5:6]=lattice[:,5:6] + multiply(dot((1 / 6),dt),(vel + dot(2.0,vel1) + dot(2.0,vel2) + vel3))
    vel=vel + multiply(dot((1 / 6),dt),(Fn + dot(2.0,F1) + dot(2.0,F2) + F3)) / m
    lattice,vel=set_boundaryconditions(tn + dt,N,lattice,vel,boundarysets,boundarytype,posfuncs,velfuncs,nargout=2)
    return lattice,vel