function[Fdamp]=evaluatedampingforces2D(t,N,lattice,dissipationparam,vel,dsthandle,dshhandle,dbehandle,structuralneighbours,shearneighbours,bendneighbours)

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

Fdamp = zeros(N,2);

%%

% ---> structural springs

for i=1:N
    for j=1:4
        neigh = structuralneighbours(i,j+1);
        if neigh~=-1
            Fdamp(i,:) = Fdamp(i,:) + dsthandle(dissipationparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum(((vel(i,:)-vel(neigh,:)).*(((lattice(neigh,5:6)-lattice(i,5:6)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))),2).*(((lattice(neigh,5:6)-lattice(i,5:6)))/(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))));
        end
    end
end

% ---> shear springs

for i=1:N
    for j=1:4
        neigh = shearneighbours(i,j+1);
        if neigh~=-1
            Fdamp(i,:) = Fdamp(i,:) + dshhandle(dissipationparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum(((vel(i,:)-vel(neigh,:)).*(((lattice(neigh,5:6)-lattice(i,5:6)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))),2).*(((lattice(neigh,5:6)-lattice(i,5:6)))/(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))));
        end
    end
end

% ---> bend springs

for i=1:N
    for j=1:4
        neigh = bendneighbours(i,j+1);
        if neigh~=-1
            Fdamp(i,:) = Fdamp(i,:) + dbehandle(dissipationparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum(((vel(i,:)-vel(neigh,:)).*(((lattice(neigh,5:6)-lattice(i,5:6)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))),2).*(((lattice(neigh,5:6)-lattice(i,5:6)))/(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))));
        end
    end
end

return