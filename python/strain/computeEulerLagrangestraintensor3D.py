# Autogenerated with SMOP 
from smop.core import *
# 

    
@function
def computeEulerLagrangestraintensor3D(C=None,*args,**kwargs):
    varargin = computeEulerLagrangestraintensor3D.varargin
    nargin = computeEulerLagrangestraintensor3D.nargin

    ##
#        Project: Fluid - structure interaction on deformable surfaces
#         Author: Luca Di Stasio
#    Institution: ETH Zrich
#                 Institute for Building Materials
# Research group: Computational Physics for Engineering Materials
#        Version: 0.1
#  Creation date: July 30th, 2014
#    Last update: July 31st, 2014
    
    #    Description: 
#          Input: 
#         Output:
    
    ##
    
    M,N=size(C,nargout=2)
    E=multiply(0.5,(C - sparse(cat([(arange(1,M)).T],[(arange(1,M)).T],[(arange(1,M)).T]),cat([ones(M,1)],[dot(4,ones(M,1))],[dot(7,ones(M,1))]),ones(dot(3,M),1),M,N)))
    return E