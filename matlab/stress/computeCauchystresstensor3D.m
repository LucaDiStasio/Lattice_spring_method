function[sigma]=computeCauchystresstensor3D(t,dt,N,lattice,oldpos,sqrtg,Vref,elasticparam,ksthandle,kshhandle,kbehandle,structuralneighbours,shearneighbours,bendneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 31st, 2014
%    Last update: July 31st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

oldFel = evaluateelasticforces3D(t-dt,N,[lattice(:,1:6) oldpos],elasticparam,ksthandle,kshhandle,kbehandle,structuralneighbours,shearneighbours,bendneighbours);

Fel    = evaluateelasticforces3D(t,N,lattice,elasticparam,ksthandle,kshhandle,kbehandle,structuralneighbours,shearneighbours,bendneighbours);

u = lattice(:,7:9) - oldpos;

sigma = 0.5.*[(oldFel(:,1)+Fel(:,1)).*u(:,1) (oldFel(:,1)+Fel(:,1)).*u(:,2) (oldFel(:,1)+Fel(:,1)).*u(:,3)...
              (oldFel(:,2)+Fel(:,2)).*u(:,1) (oldFel(:,2)+Fel(:,2)).*u(:,2) (oldFel(:,2)+Fel(:,2)).*u(:,3)...
              (oldFel(:,3)+Fel(:,3)).*u(:,1) (oldFel(:,3)+Fel(:,3)).*u(:,2) (oldFel(:,3)+Fel(:,3)).*u(:,3)]./[(sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref) (sqrtg.*Vref)];

return