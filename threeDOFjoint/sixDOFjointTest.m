out = timerfind;
delete(out);

clc; close all; clear all;

addpath('..\tensegrityObjects')
r1 = 0.15;  %radius bot tetrahedron
h1 = 0.2;
h2 = 0.1652;%0.0436;%  %radius top tetrahedron
r2 = 0.0478;%%height top
NR = 0.875; %ratio of nesting

quat = dcm2quat(getHG_Tform(0,0,0));
jointNodes = get3DOFJointNodes(h1,r1,r2,h2,NR,quat);
T= 0 ; G = 0; P = 0;

stringStiffness = [300000 * ones(3,1); 5000*ones(6,1)];
barStiffness = 300000*ones(13,1);
stringDamping = 10000*ones(9,1);
barDamping = 50000*ones(13,1);
nodalMass = 0.1*ones(8,1);
nodalMass(end) = 2;
nodalMass(1:3) = 0;
g=9.81; 
F = [zeros(8,1), zeros(8,1), -nodalMass*9.81];
F(1:3,3) = -sum(F(:,3))*ones(3,1)/3;

tspan = 0.01;

delT = 0.0005;
minQ = 0;

lims = 0.3;

bars = [1 2 3 4 4 4 5 6 7 8 8 8 8;
        2 3 1 5 6 7 6 7 5 5 6 7 4];
    
strings = [1 2 3 1 1 2 2 3 3;
           4 4 4 7 5 5 6 6 7];

%%%%%%%%%%%%%% Create objects for spine IK/dynamics, and plotting objects %%%%%%
% Pass all variables to the TensegrityStructure constructor to create the
% object which we call spine. This object contains methods for inverse
% kinematics as well as dynamics
joint = TensegrityStructure(jointNodes, strings, bars, F,stringStiffness,...
                           barStiffness,stringDamping,nodalMass,delT,[],[]);
 
theta = (0:(2*pi/20):2*pi)';
ringNodes = [r1*sin(theta), r1*cos(theta), jointNodes(1,3)*ones(size(theta)) ];
 ringStrings = [1;3];
 ringBars = [1:length(theta);
             length(theta), 1:length(theta)-1];
     
jointCommandPlot = TensegrityPlot(jointNodes, strings, bars(:,4:end), 0.005,0.001);
jointDynamicsPlot = TensegrityPlot(jointNodes, strings, bars(:,4:end), 0.005,0.001);
ringPlot = TensegrityPlot(ringNodes, ringStrings, ringBars, 0.005,0.001);

f = figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%% IK Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = subplot(1,2,1,'Parent',f,'units','normalized','outerposition',...
    [0.01 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(jointCommandPlot,ax)
updatePlot(jointCommandPlot)

generatePlot(ringPlot,ax)
updatePlot(ringPlot)
%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 1],'Style','local')
%lighting flat
lighting gouraud
colormap cool% winter
xlim([-lims lims])
ylim([-lims lims])
zlim([-0.01 1.6*lims])
title('Static IK Command');

%%%%%% Dynamics Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = subplot(1,2,2,'Parent',f,'units','normalized','outerposition',...
    [0.51 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(jointDynamicsPlot,ax)
updatePlot(jointDynamicsPlot)

generatePlot(ringPlot,ax)
updatePlot(ringPlot)

axis equal
view(3)
grid on
light('Position',[0 0 1],'Style','local')
xlim([-lims lims])
ylim([-lims lims])
zlim([-0.01 1.6*lims])
title('Dynamic Response');
% Create an object calls from the class TensegrityCallbackFunctions which is
% just a function wrapper to hold a vector of persistent values, and some update functions
% It also holds a function for creating sliders
calls = TensegrityCallbackFunctions([T, G ,P]);

%Use the function to create some slider below and attach some callback
%functions to update angle, bed axis, twist, and NR
bendLims = 1.1;
calls.makeSlider(f,0,0,@(es,ed) updateVal(calls,es.Value,1),[-bendLims bendLims],'Theta')
calls.makeSlider(f,475,0,@(es,ed) updateVal(calls,es.Value,2),[-bendLims bendLims],'Gamma')
calls.makeSlider(f,2*475,0,@(es,ed) updateVal(calls,es.Value,3),[-0.8 0.8],'Phi')

pStruct = struct('minQ',minQ,'tspan',tspan,'r1',r1,'h1',h1,'r2',r2,'h2',h2,'NR',NR);
vec1 = [T,G,P];
jointUpdate(vec1,...
    joint, jointCommandPlot,jointDynamicsPlot,pStruct);

%Create a function handle which only passes the vector of values we will be
%updating
jointUpdates = @(vec) jointUpdate(vec);
%Create a timer to update everything, 20 fps should
%look smooth prob best not to go below this
% for i = 1:1000
%    timerUpdate(calls,jointUpdates)
% end
t = timer;
t.TimerFcn = @(myTimerObj, thisEvent) timerUpdate(calls,jointUpdates);
t.Period = tspan;
t.ExecutionMode = 'fixedRate';
start(t);




       