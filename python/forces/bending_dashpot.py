# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def bending_dashpot(dissipationparam=None,t=None,points=None,lattice=None,*args,**kwargs):
    varargin = bending_dashpot.varargin
    nargin = bending_dashpot.nargin

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
    
    Q=dissipationparam[1]
    E=dissipationparam[2]
    d=E / Q
    return d