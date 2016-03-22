function[Fselfcontact]=evaluateselfcontact3D(N,lattice,boundary,Rcsq,ksc)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 28th, 2014
%    Last update: July 29th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Nb = size(boundary,1);

indicesrow = [];
indicescolumn = [];
distances = [];

index = (1:Nb)';

for row=1:Nb-1
    columns = index(row+1:end);
    numcol = length(columns);
    indicesrow = [indicesrow;row*ones(length(columns),1)];
    indicescolumn = [indicescolumn;columns];
    distances = [distances; sum(([lattice(boundary(row),7)*ones(length(columns),1) lattice(boundary(row),8)*ones(length(columns),1) lattice(boundary(row),9)*ones(length(columns),1)]-lattice(boundary(columns),7:9)).^2,2)];
end

clear val numcol columns row

[rows,cols] = find(sparse(indicesrow,indicescolumn,distances,Nb,Nb)<Rcsq);

clear indicesrow indicescolumn distances

if ~isempty(rows)
    F = sparse((1:3*N)',ones(3*N,1),zeros(3*N,1),3*N,1);
    for i=1:length(rows)
        Fc2 = ksc.*(lattice(boundary(cols(i)),7:9)-lattice(boundary(rows(i)),7:9));
        Fc1 = -Fc2;
        F(boundary(rows(i)),1)       = F(boundary(rows(i)),1)       + Fc1(1,1);
        F(N + boundary(rows(i)),1)   = F(N + boundary(rows(i)),1)   + Fc1(1,2);
        F(2*N + boundary(rows(i)),1) = F(2*N + boundary(rows(i)),1) + Fc1(1,3);
        F(boundary(cols(i)),1)       = F(boundary(cols(i)),1)       + Fc2(1,1);
        F(N + boundary(cols(i)),1)   = F(N + boundary(cols(i)),1)   + Fc2(1,2);
        F(2*N + boundary(cols(i)),1) = F(2*N + boundary(cols(i)),1) + Fc2(1,3);
    end
    clear rows cols
    Fselfcontact = sparse([(1:N)';(1:N)';(1:N)'],[ones(N,1);2*ones(N,1);3*ones(N,1)],F,N,3);
else
    clear rows cols
    Fselfcontact = sparse([(1:N)';(1:N)';(1:N)'],[ones(N,1);2*ones(N,1);3*ones(N,1)],zeros(3*N,1),N,3);
end

return