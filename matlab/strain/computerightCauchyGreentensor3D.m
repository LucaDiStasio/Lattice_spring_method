function[C]=computerightCauchyGreentensor3D(F)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 30th, 2014
%    Last update: July 31st, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

C = [F(:,1).*F(:,1)+F(:,4).*F(:,4)+F(:,7).*F(:,7) F(:,1).*F(:,2)+F(:,4).*F(:,5)+F(:,7).*F(:,8) F(:,1).*F(:,3)+F(:,4).*F(:,6)+F(:,7).*F(:,9)...
     F(:,2).*F(:,1)+F(:,5).*F(:,4)+F(:,8).*F(:,7) F(:,2).*F(:,2)+F(:,5).*F(:,5)+F(:,8).*F(:,8) F(:,2).*F(:,3)+F(:,5).*F(:,6)+F(:,8).*F(:,9)...
     F(:,3).*F(:,1)+F(:,6).*F(:,4)+F(:,9).*F(:,7) F(:,3).*F(:,2)+F(:,6).*F(:,5)+F(:,9).*F(:,8) F(:,3).*F(:,3)+F(:,6).*F(:,6)+F(:,9).*F(:,9)];

return