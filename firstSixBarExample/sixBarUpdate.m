function myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval)
% This function will perform dynamics update each timestep.

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan i

if nargin>1
    i = 0;
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
end

%%% Optional rest-length controller %%%
i = i + 1;
if i == 50 % Start after a certain time.
    newRestLengths = rand(24, 1); % Random rest lengths.
    tensStruct.simStruct.stringRestLengths = newRestLengths;
end

%%% End controller %%%

% Update nodes:
dynamicsUpdate(tensStruct, tspan);
dynamicsPlot.nodePoints = tensStruct.ySim(1:end/2,:);
updatePlot(dynamicsPlot);

drawnow  %plot it up
end

