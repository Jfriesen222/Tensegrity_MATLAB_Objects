function [nodeXYZ] = get3DOFJointNodes(h1,r1,r2,h2,NR,quats)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% r1: radius of base tetrahedron
% r2: radius of top  tetrahedron
% h1: height of base tetrahedron
% h2: height of top  tetrahedron
% NR: Nested Ratio,  the top tetrahedron point will be a certain ratio of
%     h1 nested into the base terahedron, in the set (0,1)
% quats: quaternion to rotate by

%Constant nodal offsets for nodes a,b,c,d (base tetrahedron)
nBase = [[r1*sin(0)       r1*cos(0)           h1*(1-NR)];
         [r1*sin(2*pi/3)  r1*cos(2*pi/3)      h1*(1-NR)];
         [r1*sin(-2*pi/3) r1*cos(-2*pi/3)     h1*(1-NR)]];
%Constant nodal offsets for nodes e,f,g,h (top tetrahedron)
nArm =  [[0               0              0];
         [r2*sin(pi/3)    r2*cos(pi/3)   h2];
         [r2*sin(pi)      r2*cos(pi)     h2];
         [r2*sin(5*pi/3)  r2*cos(5*pi/3) h2];
         [0               0              2*h2]];
%apply rotations to top tetrahedron
nArm = quatrotate(quats,nArm);
nodeXYZ = [nBase;
           nArm];
       nodeXYZ(:,3) = nodeXYZ(:,3) + 0.1;
end

