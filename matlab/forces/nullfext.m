function[F]=nullfext(pos,vel,t)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 29th, 2014
%    Last update: July 29th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

N = size(pos,1);

F = sparse([(1:N)';(1:N)';(1:N)'],[ones(N,1);2*ones(N,1);3*ones(N,1)],zeros(3*N,1),N,3);

return