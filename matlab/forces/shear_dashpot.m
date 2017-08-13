function[d]=shear_dashpot(dissipationparam,t,points,lattice)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: August 4th, 2014
%    Last update: August 4th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Q  = dissipationparam(1);
E  = dissipationparam(2);
nu  = dissipationparam(3);

G = 0.5*E/(1+nu);

d = G/Q;

return