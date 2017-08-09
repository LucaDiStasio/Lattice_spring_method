# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def shear_spring(elasticparam=None,t=None,points=None,lattice=None,*args,**kwargs):
    varargin = shear_spring.varargin
    nargin = shear_spring.nargin

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
    
    E=elasticparam[1]
    nu=elasticparam[2]
    G=dot(0.5,E) / (1 + nu)
    k=copy(G)
    return k