function[k]=shear_spring(elasticparam,t,points,lattice)

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

E  = elasticparam(1);
nu = elasticparam(2);

G = 0.5*E/(1+nu);

k = G;

return