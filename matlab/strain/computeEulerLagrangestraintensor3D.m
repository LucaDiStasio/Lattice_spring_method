function[E]=computeEulerLagrangestraintensor3D(C)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 30th, 2014
%    Last update: July 31st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

[M,N] = size(C);

E = 0.5.*(C-sparse([(1:M)';(1:M)';(1:M)'],[ones(M,1);4*ones(M,1);7*ones(M,1)],ones(3*M,1),M,N));

return