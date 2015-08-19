function superBallUpdate(superBall1,superBallDynamicsPlot1,superBallUKFPlot1,tspan1,ax1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%create some persistent variables for objects and structs
persistent superBall superBallDynamicsPlot superBallUKFPlot tspan allMeasureIndices lambdaErrors bars i ax

if nargin>1
    i = 0;
    superBall = superBall1;
    superBallDynamicsPlot = superBallDynamicsPlot1;
    superBallUKFPlot = superBallUKFPlot1;
    tspan = tspan1;
    allMeasureIndices = [2*ones(1,1), 3*ones(1,2), 4*ones(1,3), 5*ones(1,4), 6*ones(1,5), ...
                7*ones(1,6), 8*ones(1,7), 9*ones(1,8), 10*ones(1,9),11*ones(1,10),12*ones(1,11)...
                13*ones(1,12), 14*ones(1,12), 15*ones(1,12), 16*ones(1,12);
                1, 1:2, 1:3, 1:4, 1:5, 1:6, 1:7, 1:8,  1:9,1:10,1:11, 1:12, 1:12, 1:12, 1:12];
    bars = [superBall.barNodes];
    allMeasureIndices(:,[1 6 15 28 45 66]) = []; %eliminate bar measures
    lambdaErrors = 0.05*randn(size(allMeasureIndices,2),1); 
    ax = ax1;
    disp(lambdaErrors);
end
i = i+1;
%%%%%%%%%%%%%%%%%% update Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% update superBall nodes %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Compute Command and update dynamics $%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mod(i,2)==0)
currentWorkingLengths = unique(round(rand(100,1)*107))+1;
else
    currentWorkingLengths = unique(round(rand(0,1)*107))+1;
end
disp(length(currentWorkingLengths))
workingMeasureIndices = allMeasureIndices(:,currentWorkingLengths);
superBall.lengthMeasureIndices = [bars workingMeasureIndices];
dynamicsUpdate(superBall,tspan);
actualNodes =  superBall.ySim(1:end/2,:);
barVec = actualNodes(bars(1,:),:) - actualNodes(bars(2,:),:);
barNorm = sqrt(barVec(:,1).^2 + barVec(:,2).^2 + barVec(:,3).^2);
barAngleFromVert = acos(barVec(:,3:3:end)./barNorm);

LI =  workingMeasureIndices;
lengthNoise = randn(size(LI,2),1)*0.05;
tiltNoise = 0.01*randn(6,1);
yyPlusBase = [actualNodes; superBall.baseStationPoints];
allVectors = (yyPlusBase(LI(1,:),:) - yyPlusBase(LI(2,:),:)).^2;
z = sqrt(sum(allVectors,2));
superBall.measurementUKFInput =[barAngleFromVert+tiltNoise; 1.7*ones(6,1); (z + lengthNoise + lambdaErrors(currentWorkingLengths))];
ukfUpdate(superBall,tspan);

superBallDynamicsPlot.nodePoints = actualNodes ;
superBallUKFPlot.nodePoints = superBall.ySimUKF(1:end/2,:);
updatePlot(superBallDynamicsPlot);
updatePlot(superBallUKFPlot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lims = 1.2*1.7;
x_Avg = mean(actualNodes(:,1));
y_Avg = mean(actualNodes(:,2));
xlim(ax(1),[-lims lims]+x_Avg)
ylim(ax(1),[-lims lims]+y_Avg)
xlim(ax(2),[-lims lims]+x_Avg)
ylim(ax(2),[-lims lims]+y_Avg)
drawnow  %plot it up
end

