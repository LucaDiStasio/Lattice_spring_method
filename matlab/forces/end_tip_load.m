function[F]=end_tip_load(pos,vel,t,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: August 4th, 2014
%    Last update: August 4th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

tmax = 1; % [s]

force = 0.1*10^3/length([indicesE9;indicesC5;indicesC6]); % [N]

N = size(pos,1);

F = sparse(N,3);

if t<=tmax
    F([indicesE9;indicesC5;indicesC6],1:3) = [zeros(length([indicesE9;indicesC5;indicesC6]),1) -(force/tmax)*t*ones(length([indicesE9;indicesC5;indicesC6]),1) zeros(length([indicesE9;indicesC5;indicesC6]),1)];
else
    F([indicesE9;indicesC5;indicesC6],1:3) = [zeros(length([indicesE9;indicesC5;indicesC6]),1) -force*ones(length([indicesE9;indicesC5;indicesC6]),1) zeros(length([indicesE9;indicesC5;indicesC6]),1)];
end
return