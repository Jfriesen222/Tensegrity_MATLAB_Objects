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
        delT      %Timestep of simulation
        
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
        groundHeight
        
        
    end
    
    methods
        function obj = TensegrityStructure(nodePoints, stringNodes, barNodes, F,stringStiffness,barStiffness,stringDamping,nodalMass,delT)
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
            if (isnumeric(barNodes) && ~any(mod(barNodes(:),1)))...
                    && size(barNodes,1) == 2
                if  (max(barNodes(:))<=obj.n) && (min(barNodes(:))>0)
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
                    
                else
                    error('barNodes entries need to be in the range of 1 to n')
                end
            else
                error('barNodes should be a 2 by bb matrix of positive integers')
            end
            
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
            topNb = obj.barNodes(1,:);
            botNb = obj.barNodes(2,:);
            topNs = obj.stringNodes(1,:);
            botNs = obj.stringNodes(2,:);
            
            
            M = ones(size(repmat(nodalMass,1,3)))./repmat(nodalMass,1,3);
            indexes = 1:length(nodalMass);
            fN = indexes(nodalMass<=0);
            barNodeXYZ =  obj.nodePoints(topNb,:) - obj.nodePoints(botNb,:);
            barLengths = sum((barNodeXYZ).^2,2).^0.5;
            obj.simStruct.topNb = obj.barNodes(1,:);
            obj.simStruct = struct('M',M,'fN',fN,'stringStiffness',stringStiffness,...
                'barStiffness',barStiffness,'C',obj.C,'barRestLengths',barLengths,'stringDamping',stringDamping,...
                'topNb',topNb,'botNb',botNb,'topNs',topNs,'botNs',botNs,'F',obj.F);
            obj.delT = delT;
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
            if(nargin>2)
                obj.ySim = y0;
            end
            if(isempty(obj.ySim))
                y = [obj.nodePoints; zeros(size(obj.nodePoints))];
            else
                y = obj.ySim;
            end
            dt = obj.delT;
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
            for i=1:round(tspan/dt)                              % calculation loop
                k_1 = getAccel(yy,yDot);
                yDot1 = yDot+k_1*(0.5*dt);
                k_2  = getAccel(yy+yDot*(0.5*dt), yDot1);
                yDot2 = yDot+k_2*(0.5*dt);
                k_3 = getAccel(yy+yDot1*(0.5*dt),yDot2);
                yDot3 = yDot+k_3*dt;
                k_4 = getAccel(yy+yDot2*dt,yDot3);
                yy = yy + (1/6)*(yDot+2*yDot1+2*yDot2+yDot3)*dt;  % main equation
                yDot = yDot + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;  % main equation
            end
            obj.ySim =[yy;yDot];
            
            function [nodeXYZdoubleDot, nodeXYZdot] = getAccel(nodeXYZ,nodeXYZdot)
                memberNodeXYZ = nodeXYZ(topN,:) - nodeXYZ(botN,:);
                memberNodeXYZdot = nodeXYZdot(topN ,:) - nodeXYZdot(botN,:);
                lengths = sqrt(sum((memberNodeXYZ).^2,2));
                memberVel = sum(memberNodeXYZ.*memberNodeXYZdot,2);
                Q = stiffness.*(restLengths ./ lengths-1) - damping.*memberVel;
                Q((isString & (restLengths>lengths | Q>0))) = 0;
                FF = CC*(memberNodeXYZ.*[Q Q Q]);
                normForces = -500000*(nodeXYZ(:,3) - groundH);
                normForces(normForces<0) = 0;
                xyDot = nodeXYZdot(:,1:2);
                xyDotMag = sqrt(sum((xyDot).^2,2));
                frictionDensity =  -0.5*normForces./xyDotMag;
                frictionDensity(xyDotMag < 0.0001) = 0;
                tangentForces = xyDot.*[frictionDensity frictionDensity];
                groundForces = [tangentForces normForces];
                nodeXYZdoubleDot = (FF+groundForces).*M;
                nodeXYZdoubleDot(:,3) = nodeXYZdoubleDot(:,3) -9.81;
                nodeXYZdoubleDot(fN,:) = 0;
            end
        end
    end
end