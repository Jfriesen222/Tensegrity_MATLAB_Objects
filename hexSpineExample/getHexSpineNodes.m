function [spineNodes] = getHexSpineNodes(r,h,NR,quats)
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
zeroNodes = [0,                                     0,                    h;
             r*sin(0:pi/3:5*pi/3)',       r*cos(0:pi/3:5*pi/3)',                   zeros(6,1)];
%Translate vertically according to NR
translateVector = [0,0,h*NR];

%determine number of tetrahedrons
N = size(quats,1);
%initialize memory for speed
spineNodes = zeros(N*7,3);
%rotate first tetrahedron accordingly
spineNodes(1:7,:) = quatrotate(quats(1,:),zeroNodes);
%rotate translation vector for shifting vector alog vertiacal axis
translateVector = quatrotate(quats(1,:),translateVector);
for i=2:N
    %figure out translation to move the bottom point of the tetrahedron to
    %the origin so that we can apply rotations about the point by shifting
    %it the origin
    shift = spineNodes((7*(i-1)-6),:);
    repShift = shift(ones(1,7),:); %replicate the single row to get 4 by 3 so we can add it to our tetrahedron nodes
    %apply tilting quaternion here and translate back to shift the
    %tetrahedron back into position and then move it along the vertical
    %axis of the previous tetrahedron to space it out
    spineNodes((7*i-6):(7*i),:) = quatrotate(quats(i,:),spineNodes((7*(i-1)-6):(7*(i-1)),:) + translateVector(ones(1,7),:) - repShift )+ repShift  ;
    %detrmine the vertical axis shift for the next tetrahedron
    translateVector = quatrotate(quats(i,:),translateVector);    
end

zeroNodes = [0,                                     0,                    -h;
             r*sin(0:pi/3:5*pi/3)',       r*cos(0:pi/3:5*pi/3)',                   zeros(6,1)];
%Translate vertically according to NR
translateVector = [0,0,-h*NR];

%determine number of tetrahedrons
N = size(quats,1);
%initialize memory for speed
spineNodes2 = zeros(N*7,3);
%rotate first tetrahedron accordingly
quats(:,1) = -quats(:,1);
spineNodes2(1:7,:) = quatrotate(quats(1,:),zeroNodes);
%rotate translation vector for shifting vector alog vertiacal axis
translateVector = quatrotate(quats(1,:),translateVector);

for i=2:N
    %figure out translation to move the bottom point of the tetrahedron to
    %the origin so that we can apply rotations about the point by shifting
    %it the origin
    shift = spineNodes2((7*(i-1)-6),:);
    repShift = shift(ones(1,7),:); %replicate the single row to get 4 by 3 so we can add it to our tetrahedron nodes
    %apply tilting quaternion here and translate back to shift the
    %tetrahedron back into position and then move it along the vertical
    %axis of the previous tetrahedron to space it out
    spineNodes2((7*i-6):(7*i),:) = quatrotate(quats(i,:),spineNodes2((7*(i-1)-6):(7*(i-1)),:) + translateVector(ones(1,7),:) - repShift )+ repShift  ;
    %detrmine the vertical axis shift for the next tetrahedron
    translateVector = quatrotate(quats(i,:),translateVector);    
end

spineNodes = [spineNodes; spineNodes2];
spineNodes(N*7 + (2:7), :) = [];
spineNodes(:,3) = spineNodes(:,3)+N*h;

end