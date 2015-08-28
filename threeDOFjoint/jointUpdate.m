function jointUpdate(vec,joint1,jointCommandPlot1,jointDynamicsPlot1,pStruct1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%create some persistent variables for objects and structs
persistent joint jointCommandPlot jointDynamicsPlot pStruct i optimOpts 

if nargin>1
    joint = joint1;
    jointCommandPlot = jointCommandPlot1;
    jointDynamicsPlot = jointDynamicsPlot1;
    pStruct = pStruct1;
    getStaticTensions(joint,10);
    i = 0;
    optimOpts = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
    %optimOpts =  optimoptions('linprog','Algorithm',  'interior-point','Display','off','TolFun',1e-3);
end
i = i+1;
%%%%%%%%%%%%%%%%%% update Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set variables from persistent structure
r1 = pStruct.r1; h1 = pStruct.h1; r2 = pStruct.r2; h2 = pStruct.h2; NR = pStruct.NR; minQ = pStruct.minQ; tspan = pStruct.tspan;

% set variables from sliders below
T = vec(1); G = vec(2); P = vec(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% compute quaternions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quat = dcm2quat(getHG_Tform(T,G,P));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% update spine nodes in command plot and spine object %%%%%%%%%%%%%%%%%%%%%
jointNodes = get3DOFJointNodes(h1,r1,r2,h2,NR,quat);
joint.nodePoints = jointNodes;
jointCommandPlot.nodePoints = jointNodes;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Compute Command and update dynamics if feasible command is generated%%%%%%%

%disp(joint.simStruct.stringRestLengths);
joint.simStruct.stringRestLengths(1:3) = 0.1521;
dynamicsUpdate(joint,tspan);
jointDynamicsPlot.nodePoints = joint.ySim(1:end/2,:);
%disp(jointDynamicsPlot.nodePoints)

barVec = jointDynamicsPlot.nodePoints(8,:)-jointDynamicsPlot.nodePoints(4,:);
twistVec = -(jointDynamicsPlot.nodePoints(8,:)+jointDynamicsPlot.nodePoints(4,:))/2 + jointDynamicsPlot.nodePoints(6,:);

barVec = barVec/norm(barVec);
twistVec = twistVec/norm(twistVec);
Gdy =  asin(barVec(1)) + randn(1)*1/180*pi;
Tdy = -asin(barVec(2)/(cos(Gdy))) + randn(1)*1/180*pi;
Pdy = asin(twistVec(1)/cos(Gdy)) + randn(1)*1/180*pi;

quatdy = dcm2quat(getHG_Tform(Tdy,Gdy,Pdy));
quatError = quatmultiply(quatdy,quatconj(quat));
minAngleError = 2*acos(quatError(1));
normAxisError = quatError(2:4)*(1/sin(minAngleError/2));
restoringMoments = [300; 300; 300].*(minAngleError*normAxisError)';
restoringMoments(isnan(restoringMoments)) = 0;
cableVectors = -jointDynamicsPlot.nodePoints([7 5 5 6 6 7],:) + jointDynamicsPlot.nodePoints([1 1 2 2 3 3],:);
mountPts = jointDynamicsPlot.nodePoints([7 5 5 6 6 7],:)  - jointDynamicsPlot.nodePoints([4 4 4 4 4 4],:);
lengths = sqrt(sum(cableVectors.^2,2));
cableVectors = (cableVectors)./lengths(:,[1 1 1]);
cableMoments = cross(mountPts,cableVectors)';

A = cableMoments;
tol = norm(A)/10;
T = quadprog( A'*A + (tol)*eye(6) , -2*restoringMoments'*A ,[],[],[],[],zeros(6,1),500*(ones(6,1)./lengths),[],optimOpts) ;
disp((T)')
%q = linprog(ones(6,1),[],[],A,restoringMoments,zeros(6,1),[],[],optimOpts);
restLengths =  (lengths-T./joint.simStruct.stringStiffness(4:end));
joint.simStruct.stringRestLengths(4:end) = restLengths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(i,10) == 0
updatePlot(jointCommandPlot);    
updatePlot(jointDynamicsPlot);
drawnow %plot it up
end
end

