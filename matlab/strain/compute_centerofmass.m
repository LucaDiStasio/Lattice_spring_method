function[cm]=compute_centerofmass(m,lattice)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 3rd, 2014
%    Last update: July 3rd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

cm = sum(m.*lattice(7:9),1)./sum(m);

return