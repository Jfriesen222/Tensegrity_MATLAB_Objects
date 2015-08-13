classdef TensegrityStructure < handle
    properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% User Set Values %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nodePoints            % n by 3 matrix of node points
        
        % Below are two matrices describing string and bar connectivity which
        % are 2 by ss and 2 by bb where ss is the number of strings and bb is
        % the number of bars. each column of this matrix corresponds to a string
        % or bar and the top and bottom entries are the node numbers that the
        % string or bar spans
        stringNodes           %2 by ss matrix of node numbers for each string
        %end node, top row must be less than bottom row
        
        barNodes              %2 by bb matrix node numbers for each bar end
        %node, top row must be less than bottom row
        F                     %n by 3 matrix nodal forces
        quadProgOptions       %options for quad prog
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Added for dynamics %%%%%%%%%%%%%%%%%%
        simStruct %a structure containing most variables needed for simulation functions
        %this improves efficiency because we don't need to pass
        %the entire object to get simulation variables
        simStructUKF
        delT      %Timestep of simulation
        delTUKF   %Timestep of UKF simulation
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% Auto Generated %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Cs                    %ss by n string connectivity matrix which is auto-generated
        Cb                    %bb by n bar connectivity matrix which is auto-generated
        C                     %(ss+bb) by n connectivity matrix which is auto-generated
        n                     %scalar number of nodes
        bb                    %scalar number of bars
        ss                    %scalar number of strings
        ySim
        
        
        ySimUKF
        groundHeight
        measurementUKFInput
        P
        lengthMeasureIndices
        baseStationPoints
        
    end
    
    methods
        function obj = TensegrityStructure(nodePoints, stringNodes, barNodes, F,stringStiffness,barStiffness,stringDamping,nodalMass,delT,delTUKF,stringRestLengths)
            if(size(nodePoints,2)~=3 || ~isnumeric(nodePoints))
                error('node points should be n by 3 matrix of doubles')
            end
            obj.nodePoints = nodePoints;
            obj.n = size(nodePoints,1);
            
            %%%%%%%%%%%%%%% Check stringNodes for errors %%%%%%%%%%%%%%%%%
            if((isnumeric(stringNodes) && ~any(mod(stringNodes(:),1)))...
                    && size(stringNodes,1) == 2 )
                if  (max(stringNodes(:))<=obj.n) && (min(stringNodes(:))>0)
                    obj.ss = size(stringNodes,2);
                    for i= 1:obj.ss
                        if stringNodes(1,i) == stringNodes(2,i)
                            error('stringnodes has identical entries in a column')
                        else if stringNodes(1,i) > stringNodes(2,i)
                                stringNodes(1:2,i) = stringNodes(2:-1:1,i);
                            end
                        end
                    end
                    obj.stringNodes = stringNodes;
                    
                else
                    error('stringNodes entries need to be in the range of 1 to n')
                end
            else
                error('stringNodes should be a 2 by ss matrix of positive integers')
            end
            
            %%%%%%%%%%%%%%% Check barNodes for errors %%%%%%%%%%%%%%%%%
            
                    obj.bb = size(barNodes,2);
                    for i= 1:obj.bb
                        if barNodes(1,i) == barNodes(2,i)
                            error('barnodes has identical entries in a column')
                        else if barNodes(1,i) > barNodes(2,i)
                                barNodes(1:2,i) = barNodes(2:-1:1,i);
                            end
                        end
                    end
                    obj.barNodes = barNodes;
                   
            
            %%%%%%%%%%%%% Check for repeat bars or strings %%%%%%%%%%%%%
            B = unique([stringNodes barNodes]', 'rows');
            if size(B,1) ~= (obj.bb+obj.ss)
                error('SOme bars or strings are repeated between node sets')
            end
            
            %%%%%%%%%%%%% Build Connectivity matrices  %%%%%%%%%%%%%%%%%%
            obj.Cs = zeros(obj.ss,obj.n);
            obj.Cb = zeros(obj.bb,obj.n);
            for i=1:obj.ss
                obj.Cs(i,stringNodes(1,i)) = 1;
                obj.Cs(i,stringNodes(2,i)) = -1;
            end
            for i=1:obj.bb
                obj.Cb(i,barNodes(1,i)) = 1;
                obj.Cb(i,barNodes(2,i)) = -1;
            end
            obj.C = ([obj.Cs; obj.Cb]);
            
            %%%%%%%%%%%%% if no forces provided set forces to zero %%%%%%%
            if( isempty(F))
                obj.F = zeros(obj.n,3);
            else
                if(size(F,1)~=obj.n || size(F,2)~=3 || ~isnumeric(nodePoints))
                    error('F should be n by 3 matrix of doubles')
                end
                obj.F = F;
            end
            
            obj.quadProgOptions = optimoptions('quadprog','Algorithm',  'interior-point-convex','Display','off');
            obj.groundHeight = 0;
            %%%%%%%%%%%%%%%%%%%%Dynamics Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %these are quick vector lists of bars and strings inserted into simstruct
            %used for efficiently computing string and bar lengths etc. so
            %that we don't waste time indexing string and bar nodes
            if(isempty(obj.barNodes))
                topNb = [];
                botNb =[];
            else
            topNb = obj.barNodes(1,:);
            botNb = obj.barNodes(2,:);
            end
            topNs = obj.stringNodes(1,:);
            botNs = obj.stringNodes(2,:);
            M = ones(size(repmat(nodalMass,1,3)))./repmat(nodalMass,1,3);
            indexes = 1:length(nodalMass);
            fN = indexes(nodalMass<=0);
            barNodeXYZ =  obj.nodePoints(topNb,:) - obj.nodePoints(botNb,:);
            barLengths = sum((barNodeXYZ.*barNodeXYZ),2).^0.5;
            %obj.simStruct.topNb = obj.barNodes(1,:);
            obj.simStruct = struct('M',M,'fN',fN,'stringStiffness',stringStiffness,...
                'barStiffness',barStiffness,'C',obj.C,'barRestLengths',barLengths,'stringDamping',stringDamping,...
                'topNb',topNb,'botNb',botNb,'topNs',topNs,'botNs',botNs,'stringRestLengths',stringRestLengths);
            nUKF =1 + 12*(obj.n );
            obj.simStructUKF = struct('nUKF',nUKF,'M',repmat(M,1,nUKF),'fN',fN,'stringStiffness',repmat(stringStiffness,1,nUKF),...
                'barStiffness',repmat(barStiffness,1,nUKF),'C',sparse(obj.C),'barRestLengths',repmat(barLengths,1,nUKF),'stringDamping',repmat(stringDamping,1,nUKF),...
                'topNb',topNb,'botNb',botNb,'topNs',topNs,'botNs',botNs,'stringRestLengths',repmat(stringRestLengths,1,nUKF));
            obj.delT = delT;
            obj.delTUKF = delTUKF;
            
        end
        
        function staticTensions = getStaticTensions(obj,minForceDensity)
            A= [obj.C' *diag(obj.C*obj.nodePoints(:,1));
                obj.C' *diag(obj.C*obj.nodePoints(:,2));
                obj.C' *diag(obj.C*obj.nodePoints(:,3))];
            A_g = pinv(A);
            A_g_A=A_g*A;
            V=(eye(size(A_g_A,2))-A_g_A);
            [V,R,~] = qr(V(1:obj.ss,:));
            R =diag(R);
            V = V(:,abs(R) > 10^-12);
            Hqp = V'*V;
            fqp = V'*A_g(1:obj.ss,:)*obj.F(:);
            Aqp = -V;
            bqp = A_g(1:obj.ss,:)*obj.F(:) - minForceDensity;
            w = quadprog(Hqp,fqp,Aqp,bqp,[],[],[],[],[],obj.quadProgOptions);
            q=A_g(1:obj.ss,:)*obj.F(:) + V*w;
            lengths = sum((obj.nodePoints(obj.simStruct.topNs,:) - obj.nodePoints(obj.simStruct.botNs,:)).^2,2).^0.5;
            setStringRestLengths(obj,q,lengths)
            staticTensions = q.*lengths;
        end
        
        function setStringRestLengths(obj,q,lengths)
            obj.simStruct.stringRestLengths =  lengths.*(1-q./obj.simStruct.stringStiffness);
        end
        
        
        function lengths = getLengths(obj,nodeXYZ)
            %Get lengths of all members
            lengths = sum((nodeXYZ([obj.simStruct.topNs obj.simStruct.topNb],:) - nodeXYZ([obj.simStruct.botNs obj.simStruct.botNb],:)).^2,2).^0.5;
        end
        
        %%%%%%%%%%%%%%%%%%% Dynamics Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dynamicsUpdate(obj,tspan,y0)
            persistent lastContact
            if(nargin>2)
                obj.ySim = y0;
            end
            if(isempty(obj.ySim))
                y = sparse([obj.nodePoints; zeros(size(obj.nodePoints))]);
                lastContact = obj.nodePoints(:,1:2);
            else
                y = obj.ySim;
            end
            dt = obj.delT;
 



            %friction model constants
            Kp = 20000;
            Kd = 5000;
            muS = 0.64;
            muD = 0.54;
            kk = 1000;
            kFP = 5000;
            kFD = 0;
            %getStateDerivative(obj.simStruct);
            sim = obj.simStruct;
            groundH = obj.groundHeight;
            M = sim.M; fN = sim.fN;
            stiffness = [sim.stringStiffness; sim.barStiffness];
            CC = sim.C';
            restLengths = [sim.stringRestLengths; sim.barRestLengths];
            damping = [sim.stringDamping; zeros(obj.bb,1)];
            topN = [sim.topNs sim.topNb];
            botN = [sim.botNs sim.botNb];
            isString = [ones(1,obj.ss) zeros(1,obj.bb)]';
            yy = y(1:end/2,:);
            yDot = y((1:end/2)+end/2,:);
            for i = 1:round(tspan/dt)                              % calculation loop
                k_1 = getAccel(yy,yDot);
                yDot1 = yDot+k_1*(1/3*dt);
                k_2  = getAccel(yy+yDot*(1/3*dt), yDot1);
                yDot2 = yDot+(k_2 - (1/3)*k_1)*(dt);
                k_3 = getAccel(yy+(yDot1-1/3*yDot)*dt,yDot2);
                yDot3 = yDot+(k_1 -k_2 + k_3)*dt;
                k_4 = getAccel(yy+(yDot-yDot1+yDot2)*dt,yDot3);
                yy = yy + (dt/8)*(yDot+3*(yDot1+yDot2)+yDot3);  % main equation
                yDot = yDot + (dt/8)*(k_1+3*(k_2+k_3)+k_4);  % main equation
                lastContact(staticNotApplied,:) = yy(staticNotApplied,1:2);
            end
            obj.ySim =[yy;yDot];
            
            function nodeXYZdoubleDot = getAccel(nodeXYZ,nodeXYZdot)
                memberNodeXYZ = nodeXYZ(topN,:) - nodeXYZ(botN,:);
                memberNodeXYZdot = nodeXYZdot(topN ,:) - nodeXYZdot(botN,:);
                lengths = sqrt(sum((memberNodeXYZ).^2,2));
                memberVel = sum(memberNodeXYZ.*memberNodeXYZdot,2);
                Q = stiffness.*(restLengths ./ lengths-1) - damping.*memberVel;
                Q((isString & (restLengths>lengths | Q>0))) = 0;                
                FF = CC*(memberNodeXYZ.*Q(:,[1 1 1]));            
                %update points not in contact
                notTouching = (nodeXYZ(:,3) - groundH)>0;
                %Compute normal forces
                normForces = (groundH-nodeXYZ(:,3)).*(Kp - Kd*nodeXYZdot(:,3)); 
                normForces(notTouching) = 0; %norm forces not touching are zero
                xyDot = nodeXYZdot(:,1:2);
                %Possible static friction to apply
                staticF = kFP*(lastContact - nodeXYZ(:,1:2)) - kFD*xyDot;
                staticNotApplied = (sum((staticF).^2,2) > (muS*normForces).^2)|notTouching;
                staticF(staticNotApplied,:) = 0;
                xyDotMag = sqrt(sum((xyDot).^2,2));
                w = (1 - exp(-kk*xyDotMag))./xyDotMag;
                w(xyDotMag<1e-9) = kk;
                dynamicFmag =  - muD * normForces .*w ;
                dynamicF = dynamicFmag(:,[1 1]).* xyDot;  
                dynamicF(~staticNotApplied,:) = 0;
                tangentForces = staticF + dynamicF ;
                groundForces = [tangentForces normForces];
                nodeXYZdoubleDot = (FF+groundForces).*M;
                nodeXYZdoubleDot(:,3) = nodeXYZdoubleDot(:,3) -9.81;
                nodeXYZdoubleDot(fN,:) = 0;
            end
        end
        
        function ukfUpdate(obj,tspan,y0)
            sim = obj.simStructUKF;
            nUKF = sim.nUKF;
            %              xyI = [1:3:3*nUKF; 2:3:3*nUKF];
            %              xyI = xyI(:);
            fIndex = [1:2:nUKF*2; 2:2:nUKF*2; (2*nUKF+1):nUKF*3];
            fIndex = fIndex(:);
            Qindex = [1:nUKF; 1:nUKF; 1:nUKF];
            Qindex = Qindex(:);
            Gindex = [1:nUKF; 1:nUKF];
            Gindex = Gindex(:);
            if(nargin>2)
                obj.ySimUKF = y0;
            end
            if(isempty(obj.ySimUKF))
                obj.ySimUKF = [obj.nodePoints; zeros(size(obj.nodePoints))];
                obj.P = eye((nUKF-1)/2);
            else
                y = obj.ySimUKF;
            end
            dt = obj.delTUKF;
            
            %%%%%%%%%%%%% ukf variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            z =  obj.measurementUKFInput(:);
            x = obj.ySimUKF(:);
            L = (nUKF-1)/2;
            LI = obj.lengthMeasureIndices;
            m = size(LI,2);
            alpha=1e-3;                                 %default, tunable
            beta=2;                                     %default, tunable
            lambda= 2-L;%alpha^2*(L+ki)-L;                    %scaling factor
            c=L+lambda;                                 %scaling factor
            Ws=[lambda/c 0.5/c+zeros(1,2*L)];
            fN = sim.fN;
            Wc=Ws;
            Wc(1) = Wc(1)+(1-alpha^2+beta^2);
            c=sqrt(c);
            X=sigmas(x,obj.P,c);
            X = reshape(X,obj.n*2,[]);
            xx = reshape(x,obj.n*2,[]);
            X(fN,:) = repmat(xx(fN,:),1,nUKF);
            X(fN+obj.n,:) = 0;
            Q_noise = 0.01^2*eye(L); %process noise covariance matrix
            R_noise = 0.05^2*eye(m); %measurement noise covariance matrix
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            groundH = obj.groundHeight;
            M = sim.M;
            stiffness = [sim.stringStiffness; sim.barStiffness];
            CC = sim.C';
            restLengths = [sim.stringRestLengths; sim.barRestLengths];
            damping = [sim.stringDamping; zeros(obj.bb,nUKF)];
            topN = [sim.topNs sim.topNb];
            botN = [sim.botNs sim.botNb];
            isString = [ones(obj.ss,nUKF); zeros(obj.bb,nUKF)];
            
            
            ind1 = 1:3:3*nUKF; ind2 = ind1+1; ind3 = ind1+2;
            
            yy = X(1:end/2,:);
            yDot = X((1:end/2)+end/2,:);
            
            for i = 1:round(tspan/dt)                              % calculation loop
                k_1 = getAccels(yy,yDot);
                yDot1 = yDot+k_1*(1/3*dt);
                k_2  = getAccels(yy+yDot*(1/3*dt), yDot1);
                yDot2 = yDot+(k_2 - (1/3)*k_1)*(dt);
                k_3 = getAccels(yy+(yDot1-1/3*yDot)*dt,yDot2);
                yDot3 = yDot+(k_1 -k_2 + k_3)*dt;
                k_4 = getAccels(yy+(yDot-yDot1+yDot2)*dt,yDot3);
                yy = yy + (dt/8)*(yDot+3*(yDot1+yDot2)+yDot3);  % main equation
                yDot = yDot + (dt/8)*(k_1+3*(k_2+k_3)+k_4);  % main equation
            end
            
            %%%%%%%%%%%%%% Unscented Transformation of Process %%%%%%%%%%%%
            
            X1 =[yy;yDot]; %Forward propagated particles
            X1 = reshape(X1,obj.n*6,[]);
            x1 = X1*Ws';    %Weighted average of forward propagated particles
            X2 = X1 - x1(:,ones(1,nUKF)); %Particles with average subtracted
            P1 = X2*diag(Wc)*X2'+Q_noise; %State covariance?
            
            
            %%%%%%%%%%%%% Unscented Transformation of Measurements %%%%%%%%
            yyPlusBase = [yy; repmat(obj.baseStationPoints,1,nUKF)];
            %disp(size(yy))
            allVectors = (yyPlusBase(LI(1,:),:) - yyPlusBase(LI(2,:),:)).^2;
            Z1 = sqrt(allVectors(:,ind1) + allVectors(:,ind2) + allVectors(:,ind3));
            % this is if you have xyz coord -> Z1 = reshape(yy,m,[]); %Measurements are just x-y-z coord for now
            z1 = Z1*Ws'; %Weighted average of forward propagated measurements
            Z2 = Z1 - z1(:,ones(1,nUKF)); %Measuremnets with average subtracted
            P2 = Z2*diag(Wc)*Z2'+R_noise; %Measurement covariance
            
            P12=X2*diag(Wc)*Z2'; %Transformed cross covariance matrix
            %P2_inv = pinv(P2);
            K=P12/P2;
            x=x1+K*(z-z1);                              %state update
            obj.P = P1 -K*P12';                                %covariance update
            obj.ySimUKF = reshape(x,[],3);
            
            function nodeXYZdoubleDot = getAccels(nodeXYZs,nodeXYZdots)
                memberNodeXYZ = nodeXYZs(topN,:) - nodeXYZs(botN,:);
                memberNodeXYZdot = nodeXYZdots(topN ,:) - nodeXYZdots(botN,:);
                memNodeXYZsq = memberNodeXYZ.^2;
                memNodeXYZdotProd = memberNodeXYZdot.* memberNodeXYZ;
                lengths = sqrt(memNodeXYZsq(:,ind1) + memNodeXYZsq(:,ind2) + memNodeXYZsq(:,ind3));
                memberVel = memNodeXYZdotProd(:,ind1) + memNodeXYZdotProd(:,ind2) + memNodeXYZdotProd(:,ind3);
                Q = stiffness.*(restLengths ./ lengths-1) - damping.*memberVel; %compute force densities
                Q((isString & (restLengths>lengths | Q>0))) = 0;
                GG = memberNodeXYZ.*Q(:,Qindex);
                FF = CC*GG;
                normForces = -25000*(nodeXYZs(:,3:3:end) - groundH);
                normForces(normForces<0) = 0;
                nodeXYZdotVert = reshape(nodeXYZdots',3,[]);
                xyDot = nodeXYZdotVert(1:2,:);
                xyDotMag = sqrt(sum((xyDot).^2,1));
                xyDotMag = reshape(xyDotMag,nUKF,[]).';
                xyDot = reshape(xyDot,2*nUKF,[])';
                frictionDensity =  -0.5*normForces./xyDotMag;
                
                
                frictionDensity(xyDotMag < 0.0001) = 0;
                tangentForces = xyDot.*frictionDensity(:,Gindex);
                groundForces = [tangentForces normForces];
                groundForces = groundForces(:,fIndex);
                nodeXYZdoubleDot = (FF+groundForces).*M;
                nodeXYZdoubleDot(:,3:3:end) = nodeXYZdoubleDot(:,3:3:end) - 9.81;
                nodeXYZdoubleDot(fN,:) = 0;
            end
        end
    end
end


function X=sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*chol(P)';
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A];
end