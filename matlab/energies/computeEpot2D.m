function[Epotp,Epot]=computeEpot2D(t,N,lattice,elasticparam,ksthandle,kshhandle,kbehandle,structuralneighbours,shearneighbours,bendneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 14th, 2014
%    Last update: August 4th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Epotp = zeros(N,1);

% ---> structural springs

for i=1:N
    for j=1:4
        neigh = structuralneighbours(i,j+1);
        if neigh~=-1
            Epotp(i,:) = Epotp(i,:) + 0.5*ksthandle(elasticparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum((1-(sqrt(sum((lattice(neigh,3:4)-lattice(i,3:4)).^2,2)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2)))).*(lattice(neigh,5:6)-lattice(i,5:6)).^2,2);
        end
    end
end

% ---> shear springs

for i=1:N
    for j=1:4
        neigh = shearneighbours(i,j+1);
        if neigh~=-1
            Epotp(i,:) = Epotp(i,:) + 0.5*kshhandle(elasticparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum((1-(sqrt(sum((lattice(neigh,3:4)-lattice(i,3:4)).^2,2)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2)))).*(lattice(neigh,5:6)-lattice(i,5:6)).^2,2);
        end
    end
end

% ---> bend springs

for i=1:N
    for j=1:4
        neigh = bendneighbours(i,j+1);
        if neigh~=-1
            Epotp(i,:) = Epotp(i,:) + 0.5*kbehandle(elasticparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum((1-(sqrt(sum((lattice(neigh,3:4)-lattice(i,3:4)).^2,2)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2)))).*(lattice(neigh,5:6)-lattice(i,5:6)).^2,2);
        end
    end
end

Epot = sum(Epotp);

return