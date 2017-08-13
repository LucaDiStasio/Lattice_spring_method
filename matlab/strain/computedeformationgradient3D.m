function[F]=computedeformationgradient3D(N,lattice,deltaq,firstdevneighbours)

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

undeformedcovariantbase = computecovariantbase3D(N,deltaq,[lattice(:,1:6) lattice(:,4:6)],firstdevneighbours);
[~,~,undeformedsqrtg] = computemetriccoefficients3D(undeformedcovariantbase);
undeformedcontravariantbase = computecontravariantbase3D(undeformedcovariantbase,undeformedsqrtg);

clear undeformedcovariantbase undeformedsqrtg

deformedcovariantbase = computecovariantbase3D(N,deltaq,lattice,firstdevneighbours);

G1 = undeformedcontravariantbase(:,1:3);
G2 = undeformedcontravariantbase(:,4:6);
G3 = undeformedcontravariantbase(:,7:9);

g1 = deformedcovariantbase(:,1:3);
g2 = deformedcovariantbase(:,4:6);
g3 = deformedcovariantbase(:,7:9);

clear undeformedcontravariantbase deformedcovariantbase

F = [g1(:,1).*G1(:,1) g1(:,1).*G1(:,2) g1(:,1).*G1(:,3)...
     g1(:,2).*G1(:,1) g1(:,2).*G1(:,2) g1(:,2).*G1(:,3)...
     g1(:,3).*G1(:,1) g1(:,3).*G1(:,2) g1(:,3).*G1(:,3)] + ...
    [g2(:,1).*G2(:,1) g2(:,1).*G2(:,2) g2(:,1).*G2(:,3)...
     g2(:,2).*G2(:,1) g2(:,2).*G2(:,2) g2(:,2).*G2(:,3)...
     g2(:,3).*G2(:,1) g2(:,3).*G2(:,2) g2(:,3).*G2(:,3)] + ...
    [g3(:,1).*G3(:,1) g3(:,1).*G3(:,2) g3(:,1).*G3(:,3)...
     g3(:,2).*G3(:,1) g3(:,2).*G3(:,2) g3(:,2).*G3(:,3)...
     g3(:,3).*G3(:,1) g3(:,3).*G3(:,2) g3(:,3).*G3(:,3)];

return