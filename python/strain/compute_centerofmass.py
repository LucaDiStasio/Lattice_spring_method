# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def compute_centerofmass(m=None,lattice=None,*args,**kwargs):
    varargin = compute_centerofmass.varargin
    nargin = compute_centerofmass.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: July 3rd, 2014
#    Last update: July 3rd, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    cm=sum(multiply(m,lattice[7:9]),1) / sum(m)
    return cm