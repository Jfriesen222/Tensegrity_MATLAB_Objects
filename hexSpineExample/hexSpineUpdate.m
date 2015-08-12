function hexSpineUpdate(vec,spine1,spineCommandPlot1,spineDynamicsPlot1,pStruct1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%create some persistent variables for objects and structs
persistent spine spineCommandPlot spineDynamicsPlot pStruct axisRot


if nargin>1
    spine = spine1;
    axisRot = 0;
    spineCommandPlot = spineCommandPlot1;
    spineDynamicsPlot = spineDynamicsPlot1;
    pStruct = pStruct1;
end

%%%%%%%%%%%%%%%%%% update Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set variables from persistent structure
r = pStruct.r; h = pStruct.h; N = pStruct.N; minQ = pStruct.minQ; tspan = pStruct.tspan;
% set variables from sliders below
angle = vec(1); axisRot = axisRot+vec(2)/10; twist = vec(3); NR = vec(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% compute quaternions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quatBend = [cos(angle), sin(angle)*[sin(axisRot) cos(axisRot) 0]];
quatTwist= [cos(twist) 0 0 sin(twist)];
quatt = quatmultiply(quatBend,quatTwist);
quats = [1 0 0 0;
         repmat(quatt,N-1,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% update spine nodes in command plot and spine object %%%%%%%%%%%%%%%%%%%%%
spineNodes = getHexSpineNodes(r,h,NR,quats);
spine.nodePoints = spineNodes;
spineCommandPlot.nodePoints = spineNodes;
updatePlot(spineCommandPlot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%Compute all reaction forces at the base of the spine%%%%%%%%%%%%%%%%%%%%%%%
F = spine.F; %retrieve F from object
moments = cross(spineNodes(8:end,:),F(8:end,:)); %compute all moments
sumMoments = sum(moments,1); %sum the moments
%do a bit of math to counteract moments, done in matrix otation for brevity
F(2:2:7,3)= ([-0.5, 0; -0.25, 0.5; -0.25, 0.5]./spineNodes(2:2:7,[2 1]))*sumMoments(1:2)'; 
F(isnan(F))=0;%some times we devide zero by zero so we set inf to 0
F(1,:) = -sum(F(2:end,:),1); %set the base node to make sum F in x y and z zero
spine.F = F;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Compute Command and update dynamics if feasible command is generated%%%%%%%
tensions = getStaticTensions(spine,minQ);
if(any(tensions<0))
    disp('IK could not find solution')
end
dynamicsUpdate(spine,tspan);
spineDynamicsPlot.nodePoints = spine.ySim(1:end/2,:);
updatePlot(spineDynamicsPlot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lims = 4/30*3.5*N/5;
x_Avg = mean(spine.ySim(1:end/2,1));
y_Avg = mean(spine.ySim(1:end/2,2));
xlim([-lims lims]+x_Avg)
ylim([-lims lims]+y_Avg)

drawnow %plot it up
end

