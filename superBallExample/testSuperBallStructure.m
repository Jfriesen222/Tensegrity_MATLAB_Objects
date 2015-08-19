
clc 
clear all 
close all

barSpacing = 0.375;
barLength = 1.7;
lims = 1.2*barLength;

tspan =0.05;          % time between plot updates in seconds
delT = 0.001;         % timestep for dynamic sim in seconds
delTUKF  = 0.005;
K = 998;              %outer rim string stiffness in Newtons/meter
nodalMass = 1.625*ones(12,1);
c = 40;             % damping constant, too lazy to figure out units.
F = zeros(12,3);
stringStiffness = K*ones(24,1);
barStiffness = 100000*ones(6,1);
stringDamping = c*ones(24,1);  %string damping vector

options = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
addpath('..\tensegrityObjects')

baseStationPoints = [0+0.96/2     ,   0-1.15/2      ,  1.63;
                         -1.362+0.96/2  ,   0-1.15/2      ,  1.6606 ;  
                         -2.4712+0.96/2  ,  1.1885-1.15/2 ,  1.9514;  
                          0.2882+0.96/2  ,  2.4010-1.15/2  ,  1.8013;  
                         -1.0626+0.96/2   , 2.4519-1.15/2 ,  1.7435 ];  
                     

                     
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
     
HH  = makehgtform('axisrotate',[1 1 0],0.3);
     nodes = (HH(1:3,1:3)*nodes')';
nodes(:,3) = nodes(:,3) +barLength*0.5+10 ;    

bars = [1:2:11; 
        2:2:12];
strings = [1  1   1  1  2  2  2  2  3  3  3  3  4  4  4  4  5  5  6  6  7  7  8  8;
           7  8  10 12  5  6 10 12  7  8  9 11  5  6  9 11 11 12  9 10 11 12  9 10];
stringRestLength = 0.9*ones(24,1)*norm(nodes(1,:)-nodes(7,:));


lengthMeasureIndices = [2*ones(1,1), 3*ones(1,2), 4*ones(1,3), 5*ones(1,4), 6*ones(1,5), ...
 7*ones(1,6), 8*ones(1,7), 9*ones(1,8), 10*ones(1,9),11*ones(1,10),12*ones(1,11)...
 13*ones(1,12), 14*ones(1,12), 15*ones(1,12), 16*ones(1,12);
 1, 1:2, 1:3, 1:4, 1:5, 1:6, 1:7, 1:8,  1:9,1:10,1:11, 1:12, 1:12, 1:12, 1:12]';
lengthMeasureIndices([1 6 15 28 45 66],:) = []; %eliminate bar measures
superBallCommandPlot = TensegrityPlot(nodes, strings, bars, 0.025,0.005);
N = 5;
% stringsMult = strings;
% barsMult = bars;
% for i = 1:N-1
%     stringsMult = [stringsMult strings+12*(i)];
%     barsMult = [barsMult bars+12*(i)];
% end
superBallDynamicsPlot = TensegrityPlot(nodes, strings, bars, 0.025,0.005);
superBall = TensegrityStructure(nodes, strings, bars, F, stringStiffness,...
    barStiffness, stringDamping, nodalMass, delT,delTUKF,stringRestLength);
superBall.baseStationPoints = baseStationPoints;
f = figure('units','normalized','outerposition',[0 0 1 1]);

%%%%%%%% IK Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = subplot(1,2,1,'Parent',f,'units','normalized','outerposition',...
    [0.01 0.1 0.48 0.9]);
% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(superBallCommandPlot,ax1)
updatePlot(superBallCommandPlot);
%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 10],'Style','local')
%lighting flat
lighting flat
colormap([0.8 0.8 1; 0 1 1])
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])
title('Simulated Input');


%%%%%% Dynamics Subplot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = subplot(1,2,2,'Parent',f,'units','normalized','outerposition',...
    [0.51 0.1 0.48 0.9]);

% use a method within TensegrityPlot class to generate a plot of the
% structure
generatePlot(superBallDynamicsPlot,ax2);
updatePlot(superBallDynamicsPlot);

%settings to make it pretty
axis equal
view(3)
grid on
light('Position',[0 0 10]);%,'Style','local')
lighting flat
xlim([-lims lims])
ylim([-lims lims])
zlim(1.6*[-0.01 lims])
title('UKF Output');
superBallUpdate(superBall,superBallCommandPlot,superBallDynamicsPlot,tspan,[ax1 ax2]);
%hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'});

for i = 1:200
    superBallUpdate
  %  MM(i) = getframe(f);
end
%filename = 'quickAnimation.avi';
%writerObj = VideoWriter(filename);
%writerObj.FrameRate = 20;
%open(writerObj);
%writeVideo(writerObj,MM);
%close(writerObj);
% 
% t = timer;
% t.TimerFcn = @(myTimerObj, thisEvent) superBallUpdate;
% t.Period = tspan;
% t.ExecutionMode = 'fixedRate';
% start(t);

% % 


