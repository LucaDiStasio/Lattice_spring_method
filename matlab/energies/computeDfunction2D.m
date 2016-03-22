function[Dfuncp,Dfunc]=computeDfunction2D(t,N,lattice,dissipationparam,vel,dsthandle,dshhandle,dbehandle,structuralneighbours,shearneighbours,bendneighbours)

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

Dfuncp = zeros(N,1);

% ---> structural springs

for i=1:N
    for j=1:4
        neigh = structuralneighbours(i,j+1);
        if neigh~=-1
            Dfuncp(i,:) = Dfuncp(i,:) + 0.5*dsthandle(dissipationparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum((sum(((vel(i,:)-vel(neigh,:)).*(((lattice(neigh,5:6)-lattice(i,5:6)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))),2).*(((lattice(neigh,5:6)-lattice(i,5:6)))/(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))).^2,2);
        end
    end
end

% ---> shear springs

for i=1:N
    for j=1:4
        neigh = shearneighbours(i,j+1);
        if neigh~=-1
            Dfuncp(i,:) = Dfuncp(i,:) + 0.5*dshhandle(dissipationparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum((sum(((vel(i,:)-vel(neigh,:)).*(((lattice(neigh,5:6)-lattice(i,5:6)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))),2).*(((lattice(neigh,5:6)-lattice(i,5:6)))/(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))).^2,2);
        end
    end
end

% ---> bend springs

for i=1:N
    for j=1:4
        neigh = bendneighbours(i,j+1);
        if neigh~=-1
            Dfuncp(i,:) = Dfuncp(i,:) + 0.5*dbehandle(dissipationparam,t,0.5*(lattice(neigh,5:6)+lattice(i,5:6)),lattice).*sum((sum(((vel(i,:)-vel(neigh,:)).*(((lattice(neigh,5:6)-lattice(i,5:6)))./(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))),2).*(((lattice(neigh,5:6)-lattice(i,5:6)))/(sqrt(sum((lattice(neigh,5:6)-lattice(i,5:6)).^2,2))))).^2,2);
        end
    end
end

Dfunc = sum(Dfuncp);

return
