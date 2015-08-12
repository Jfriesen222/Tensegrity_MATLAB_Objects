function [spineNodes] = getSpineNodes(r,h,NR,quats)
% % % Create an array of nodes for a tensegrity spinal column
% % 
% % %inputs:
% % % r - radius of a circle that is inscribed by th top triangle of the tetrahedra 
% % % h - height of the tetrahedra
% % % NR - nested ratio of each tetrahedron in the range 0 to 1  
% % %      h_offset = (NR*h)
% % % quats- an N by 4 matrix whose rows are quaternions and N is the number of
% % %        tetrahedra
% % 
% % %outputs:
% % %spineNodes - an N by 3 matrix whose rows re the cartesian coordinates of
% % %             each tetrahedron in space

%Create a matrix of a standar tetrahedron
zeroNodes = [0,              0,                    0;
    r*sin(0),       r*cos(0),             h;
    r*sin(2*pi/3),  r*cos(2*pi/3),        h;
    r*sin(-2*pi/3), r*cos(-2*pi/3),       h];
%Translate vertically according to NR
translateVector = [0,0,h*NR];

%need z axis because each succesive tetrahedron is rotated
% by 60 degrees for this topology
spinAxis = [0,0,1];
%determine number of tetrahedrons
N = size(quats,1);
%initialize memory for speed
spineNodes = zeros(N*4,3);
%rotate first tetrahedron accordingly
spineNodes(1:4,:) = quatrotate(quats(1,:),zeroNodes);
%rotate translation vector for shifting vector alog vertiacal axis
translateVector = quatrotate(quats(1,:),translateVector);
for i=2:N
    %figure out translation to move the bottom point of the tetrahedron to
    %the origin so that we can apply rotations about the point by shifting
    %it the origin
    shift = spineNodes((4*(i-1)-3),:);
    repShift = repmat(shift,4,1); %replicate the single row to get 4 by 3 so we can add it to our tetrahedron nodes
    
    %detrmine what the local vertical axis is to spin each tetrahedron 60 degrees as described above
    spinAxis = quatrotate(quats(i-1,:),spinAxis);
    
    %apply 60 degree spin here and also shift bottom point to the origin
    spineNodes((4*i-3):(4*i),:) = quatrotate([cos(-pi/6), sin(-pi/6)*spinAxis],...
                                             spineNodes((4*(i-1)-3):(4*(i-1)),:) - repShift );    
    
    %apply tilting quaternion here and translate back to shift the
    %tetrahedron back into position and then move it along the vertical
    %axis of the previous tetrahedron to space it out
    spineNodes((4*i-3):(4*i),:) = quatrotate(quats(i,:),spineNodes((4*(i)-3):(4*(i)),:) )+ repShift  + repmat(translateVector,4,1);
    
    %detrmine the vertical axis shift for the next tetrahedron
    translateVector = quatrotate(quats(i,:),translateVector);    
end
end