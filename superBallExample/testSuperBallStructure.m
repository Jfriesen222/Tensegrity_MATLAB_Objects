
clc 
clear variables 
close all
barSpacing = 0.375;
barLength = 1.7;
lims = 1.2*barLength;
addpath('..\tensegrityObjects')
tspan =0.045;          % time between plot updates in seconds
minQ = 100;
delT = 0.001;         % timestep for dynamic sim in seconds
K = 998;             %outer rim string stiffness in Newtons/meter
nodalMass = 1.625*ones(12,1);
c = 50;             % damping constant, too lazy to figure out units.
F = zeros(12,3);
stringStiffness = K*ones(24,1);
barStiffness = 100000*ones(6,1);
stringDamping = c*ones(24,1);  %string damping vector

options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');


nodes = [-barSpacing     barLength*0.5  0;
         -barSpacing    -barLength*0.5  0;
          barSpacing     barLength*0.5  0;
          barSpacing    -barLength*0.5  0;
          0             -barSpacing     barLength*0.5;
          0             -barSpacing    -barLength*0.5;
          0              barSpacing     barLength*0.5;
          0              barSpacing    -barLength*0.5;        
          barLength*0.5  0             -barSpacing;
         -barLength*0.5  0             -barSpacing;
          barLength*0.5  0              barSpacing;
         -barLength*0.5  0              barSpacing];
     
HH  = makehgtform('axisrotate',[1 1 0],0.1);
     nodes = (HH(1:3,1:3)*nodes')';
nodes(:,3) = nodes(:,3) +barLength*0.5 + 8;    

bars = [1:2:11; 2:2:12];
strings = [1  1   1  1  2  2  2  2  3  3  3  3  4  4  4  4  5  5  6  6  7  7  8  8;
           7  8  10 12  5  6 10 12  7  8  9 11  5  6  9 11 11 12  9 10 11 12  9 10];
stringRestLength = 0.95*ones(24,1)*norm(nodes(1,:)-nodes(7,:));
%stringRestLength(1) = stringRestLength(1)*0.9;
         
superBallCommandPlot = TensegrityPlot(nodes, strings, bars, 0.025,0.005);
superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, 0.025,0.005);
superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT,[],[]);
superBall.simStruct.stringRestLengths = stringRestLength;

f = figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%% IK Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = subplot(1,2,1,'Parent',f,'units','normalized','outerposition',...
    [0.01 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(superBallCommandPlot,ax)
updatePlot(superBallCommandPlot);
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
zlim(1.6*[-0.01 lims])
title('Static IK Command');


%%%%%% Dynamics Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = subplot(1,2,2,'Parent',f,'units','normalized','outerposition',...
    [0.51 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(superBallDynamicsPlot,ax);
updatePlot(superBallDynamicsPlot);

x = lims*[-1 1 1 -1];
y = lims*[-1 -1 1 1];
patch(x,y,zeros(size(x)))
%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 1],'Style','local')
lighting flat
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])
title('Dynamics Simulation');
pStruct = struct('minQ',minQ,'tspan',tspan);
superBallUpdate(superBall,superBallCommandPlot,superBallDynamicsPlot,pStruct);
% 
% for i = 1:200
%     superBallUpdate
% end

t = timer;
t.TimerFcn = @(myTimerObj, thisEvent) superBallUpdate;
t.Period = tspan;
t.ExecutionMode = 'fixedRate';
start(t);




