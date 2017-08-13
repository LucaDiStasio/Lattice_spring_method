function[Fextcontact]=evaluateexternalcontact3D(t,N,lattice,boundary,boundhandle,gradienthandle,kec)

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
%    Description: F(r,t) <= 0 solid obstacle
%          Input: 
%         Output: 

%%

boundindices = find(boundhandle(lattice(boundary,7:9),t)<=0);

if ~isempty(boundindices)
    F = sparse((1:3*N)',ones(3*N,1),zeros(3*N,1),3*N,1);
    for i=1:length(boundindices)
        Fc = kec.*sqrt(abs(boundhandle(lattice(boundary(boundindices(i)),7:9),t))).*gradienthandle(lattice(boundary(boundindices(i)),7:9),t)./sqrt(sum((gradienthandle(lattice(boundary(boundindices(i)),7:9),t)).^2,2));
        F(boundary(boundindices(i)),1)       = F(boundary(boundindices(i)),1)       + Fc(1,1);
        F(N + boundary(boundindices(i)),1)   = F(N + boundary(boundindices(i)),1)   + Fc(1,2);
        F(2*N + boundary(boundindices(i)),1) = F(2*N + boundary(boundindices(i)),1) + Fc(1,3);
    end
    clear boundindices
    Fextcontact = sparse([(1:N)';(1:N)';(1:N)'],[ones(N,1);2*ones(N,1);3*ones(N,1)],F,N,3);
else
    clear boundindices
    Fextcontact = sparse([(1:N)';(1:N)';(1:N)'],[ones(N,1);2*ones(N,1);3*ones(N,1)],zeros(3*N,1),N,3);
end

return