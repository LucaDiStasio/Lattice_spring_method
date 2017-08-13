
%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 7th, 2014
%    Last update: July 7th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

clear all
close all
clc

r0 = 3;
%r = r0 + (0:0.01:15)';
r = r0 + [(2:0.01:5)';(5:-0.01:2.5)';(2.5:0.01:6)';(6:-0.01:3)';(3:0.01:9)';(9:-0.01:0)';(0:0.01:11)';(11:-0.01:8.5)';(8.5:0.01:15)';(15:-0.01:0)';(0:0.01:20)'];
v = [r(2:end,1)-r(1:end-1,1);1];

delta = r - r0;
deltadot = 2*v.*(delta);
F = zeros(length(delta),1);

deltay = 3;
deltacr = 7;
deltad = 12;
deltah = 1;

Ke0 = 5;
Kp = 1.5;
Ks = (Ke0*deltay+Kp*(deltacr-deltay))/(deltad - deltacr);
Kheal = Ke0/(deltad-deltah);

delta0 = 0;
Ke = Ke0;
broken = 0;
Kehist = [Ke0];

figure
grid on
hold on
xlabel('\Delta [mm]')
ylabel('F [N]')
title('Elasto-plasticity with damage')
for i=1:length(delta)
    if ~broken
        if delta(i)>=delta0 && delta(i)<delta0+deltay
            F(i) = Ke*(delta(i)-delta0);
            plot(delta(i),F(i),'*b','LineWidth',2)
            hold on
        elseif delta(i)>=delta0+deltay && delta(i)<deltacr && deltadot(i)>=0
            F(i) = Ke0*deltay + Kp*(delta(i)-deltay);
            plot(delta(i),F(i),'*r','LineWidth',2)
            hold on
        elseif delta(i)>=delta0+deltay && delta(i)<deltacr && deltadot(i)<0
            delta0 = delta(i) - deltay;
            Ke = Ke0 + Kp*(delta0/deltay);
            Kehist = [Kehist;Ke];
            F(i) = Ke*(delta(i)-delta0);
            plot(delta(i),F(i),'or','LineWidth',2)
            hold on
        elseif delta(i)>=deltacr && delta(i)<deltad && delta(i)-delta0>=deltay && deltadot(i)>=0
            F(i) = Ke0*deltay + Kp*(deltacr-deltay) - Ks*(delta(i)-deltacr);
            plot(delta(i),F(i),'*g','LineWidth',2)
            hold on
        elseif delta(i)>=deltacr && delta(i)<deltad && deltadot(i)<0
            delta0 = delta(i) - deltay;
            Ke = Ke0 + Kp*(deltacr/deltay - 1) - Ks*(delta0/deltay-(deltacr/deltay-1));
            Kehist = [Kehist;Ke];
            F(i) = Ke*(delta(i)-delta0);
            plot(delta(i),F(i),'og','LineWidth',2)
            hold on
        elseif delta(i) == deltad
            delta0 = 0;
            Ke = 0;
            Kehist = [Kehist;Ke];
            F(i) = 0;
            broken = 1;
            plot(delta(i),F(i),'*k','LineWidth',2)
            hold on
        end
    else
        if delta(i) > deltad
            F(i) = 0;
            plot(delta(i),F(i),'*k','LineWidth',2)
            hold on
        elseif delta(i)<=deltad && delta(i)>deltah && deltadot(i)<0
            Ke = Kheal*(deltad-delta(i));
            Kehist = [Kehist;Ke];
            F(i) = Ke0*(deltah-delta0)*(Ke*(deltad-delta(i))/(Kheal*(deltad-deltah)^2));
            plot(delta(i),F(i),'*c','LineWidth',2)
            hold on
        elseif delta(i)==deltah && deltadot(i)<0
            broken = 0;
            Ke = Ke0;
            Kehist = [Kehist;Ke];
            F(i) = Ke*(delta(i)-delta0);
            plot(delta(i),F(i),'*c','LineWidth',2)
            hold on
        end
    end
end



