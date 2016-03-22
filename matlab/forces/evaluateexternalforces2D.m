function[F]=evaluateexternalforces2D(t,lattice,vel,fexthandle)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 14th, 2014
%    Last update: July 14th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

F = fexthandle(lattice(:,5:6),vel,t);

return