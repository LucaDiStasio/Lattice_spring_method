function[Ekinpdof,Ekinp,Ekin]=computeEkin(m,vel)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 2nd, 2014
%    Last update: July 2nd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Ekinpdof = 0.5.*vel.^2./[m m m];

Ekinp = 0.5.*sum(vel.^2,2)./m;

Ekin = 0.5.*sum(sum(vel.^2,2))./m;

return