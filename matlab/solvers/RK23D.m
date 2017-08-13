function[lattice,vel]=RK23D(elasticparam,dissipationparam,tn,dt,N,Nx,Ny,lattice,vel,m,ksthandle,kshhandle,kbehandle,dsthandle,dshhandle,dbehandle,fexthandle,boundarysets,boundarytype,posfuncs,velfuncs,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 2nd, 2014
%    Last update: August 4th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

Fel = evaluateelasticforces3D(tn,N,Nx,Ny,lattice,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

Fdamp = evaluatedampingforces3D(tn,N,Nx,Ny,lattice,dissipationparam,vel,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

Fext = evaluateexternalforces3D(tn,lattice,vel,fexthandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

x1 =  lattice(:,7:9) + 0.5*dt.*vel;
vel1 = vel + 0.5*dt.*(Fel+Fdamp+Fext)./m;

[lattice1,vel1] = set_boundaryconditions(tn+0.5*dt,N,[lattice(:,1:6) x1],vel1,boundarysets,boundarytype,posfuncs,velfuncs);

Fel = evaluateelasticforces3D(tn+0.5*dt,N,Nx,Ny,lattice1,elasticparam,ksthandle,kshhandle,kbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

Fdamp = evaluatedampingforces3D(tn+0.5*dt,N,Nx,Ny,lattice1,dissipationparam,vel1,dsthandle,dshhandle,dbehandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

Fext = evaluateexternalforces3D(tn+0.5*dt,lattice1,vel1,fexthandle,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8);

lattice(:,7:9) =  lattice(:,7:9) + dt.*vel1;
vel = vel + dt.*(Fel+Fdamp+Fext)./m;

[lattice,vel] = set_boundaryconditions(tn+dt,N,lattice,vel,boundarysets,boundarytype,posfuncs,velfuncs);

return