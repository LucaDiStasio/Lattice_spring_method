# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def end_tip_load(pos=None,vel=None,t=None,indicesbulk=None,indicesinternalbulk=None,indicesF1=None,indicesF2=None,indicesF3=None,indicesF4=None,indicesF5=None,indicesF6=None,indicesinternalF1=None,indicesinternalF2=None,indicesinternalF3=None,indicesinternalF4=None,indicesinternalF5=None,indicesinternalF6=None,indicesE1=None,indicesE2=None,indicesE3=None,indicesE4=None,indicesE5=None,indicesE6=None,indicesE7=None,indicesE8=None,indicesE9=None,indicesE10=None,indicesE11=None,indicesE12=None,indicesinternalE1=None,indicesinternalE2=None,indicesinternalE3=None,indicesinternalE4=None,indicesinternalE5=None,indicesinternalE6=None,indicesinternalE7=None,indicesinternalE8=None,indicesinternalE9=None,indicesinternalE10=None,indicesinternalE11=None,indicesinternalE12=None,indicesC1=None,indicesC2=None,indicesC3=None,indicesC4=None,indicesC5=None,indicesC6=None,indicesC7=None,indicesC8=None,indicesinternalC1=None,indicesinternalC2=None,indicesinternalC3=None,indicesinternalC4=None,indicesinternalC5=None,indicesinternalC6=None,indicesinternalC7=None,indicesinternalC8=None,*args,**kwargs):
    varargin = end_tip_load.varargin
    nargin = end_tip_load.nargin

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
    
    tmax=1
    
    force=dot(0.1,10 ** 3) / length(cat([indicesE9],[indicesC5],[indicesC6]))
    
    N=size(pos,1)
    F=sparse(N,3)
    if t <= tmax:
        F[cat([indicesE9],[indicesC5],[indicesC6]),1:3]=cat(zeros(length(cat([indicesE9],[indicesC5],[indicesC6])),1),dot(dot(- (force / tmax),t),ones(length(cat([indicesE9],[indicesC5],[indicesC6])),1)),zeros(length(cat([indicesE9],[indicesC5],[indicesC6])),1))
    else:
        F[cat([indicesE9],[indicesC5],[indicesC6]),1:3]=cat(zeros(length(cat([indicesE9],[indicesC5],[indicesC6])),1),dot(- force,ones(length(cat([indicesE9],[indicesC5],[indicesC6])),1)),zeros(length(cat([indicesE9],[indicesC5],[indicesC6])),1))
    
    return F