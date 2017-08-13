function[Fel]=evaluateelasticforces3D(t,N,lattice,elasticparam,ksthandle,kshhandle,kbehandle,structuralneighbours,shearneighbours,bendneighbours)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: June 26th, 2014
%    Last update: August 4th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Fel = zeros(N,3);

%%

% ---> structural springs

for i=1:N
    for j=1:6
        neigh = structuralneighbours(i,j+1);
        if neigh~=-1
            Fel(i,:) = Fel(i,:) + ksthandle(elasticparam,t,0.5*(lattice(neigh,7:9)+lattice(i,7:9)),lattice).*(1-(sqrt(sum((lattice(neigh,4:6)-lattice(i,4:6)).^2,2)))./(sqrt(sum((lattice(neigh,7:9)-lattice(i,7:9)).^2,2)))).*(lattice(neigh,7:9)-lattice(i,7:9));
        end
    end
end

% ---> shear springs

for i=1:N
    for j=1:20
        neigh = shearneighbours(i,j+1);
        if neigh~=-1
            Fel(i,:) = Fel(i,:) + kshhandle(elasticparam,t,0.5*(lattice(neigh,7:9)+lattice(i,7:9)),lattice).*(1-(sqrt(sum((lattice(neigh,4:6)-lattice(i,4:6)).^2,2)))./(sqrt(sum((lattice(neigh,7:9)-lattice(i,7:9)).^2,2)))).*(lattice(neigh,7:9)-lattice(i,7:9));
        end
    end
end

% ---> bend springs

for i=1:N
    for j=1:6
        neigh = bendneighbours(i,j+1);
        if neigh~=-1
            Fel(i,:) = Fel(i,:) + kbehandle(elasticparam,t,0.5*(lattice(neigh,7:9)+lattice(i,7:9)),lattice).*(1-(sqrt(sum((lattice(neigh,4:6)-lattice(i,4:6)).^2,2)))./(sqrt(sum((lattice(neigh,7:9)-lattice(i,7:9)).^2,2)))).*(lattice(neigh,7:9)-lattice(i,7:9));
        end
    end
end

return