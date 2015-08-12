classdef TensegrityPlot < handle
    properties
        nodePoints            % n by 3 matrix of node points
        stringNodes           %2 by ss matrix of node numbers for each string
                              %end node, top row must be less than bottom row
        barNodes              %2 by bb matrix node numbers for each bar end
                              %node, top row must be less than bottom row
        n                     %scalar number of nodes
        bb                    %scalar number of bars
        ss                    %scalar number of strings
        barRad                %bar radius for plotting
        stringRad             %string radius for plotting
        sphereTForm           %transform object for spheres which sit at nodes in plot
        barTForm              %transform object for bar cylinders
        stringTForm           %transform object for string cylinders
    end
    methods
        function obj = TensegrityPlot(nodePoints, stringNodes, barNodes, barRad, stringRad)
            if(size(nodePoints,2)~=3 || ~isnumeric(nodePoints))
                error('node points should be n by 3 matrix of doubles')
            end
            obj.nodePoints = nodePoints;
            obj.n = size(nodePoints,1);
            if(isscalar(barRad) && barRad>0)
                obj.barRad = barRad;
            else
                error('barRad requires positive scalar double value')
            end
            
            if(isscalar(stringRad) && stringRad>0)
                obj.stringRad = stringRad;
            else
                error('stringRad requires positive scalar double value')
            end
            
            
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
        end
        
        %%%%%%%%%%%%%%%%%%%% Plotting functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generatePlot(obj,ax)
            [xx,yy,zz] = sphere(); xx = xx*obj.barRad; yy = yy*obj.barRad; zz = zz*obj.barRad;
            obj.sphereTForm = gobjects(obj.n,1);
            hold(ax,'on');
            for i = 1:obj.n
                sphereHandle = surf(ax,xx,yy,zz,'LineStyle', 'none');
                obj.sphereTForm(i) = hgtransform('Parent',ax);
                set(sphereHandle,'Parent',obj.sphereTForm(i))
                set(obj.sphereTForm(i),'matrix',makehgtform('translate',obj.nodePoints(i,:)));
            end
            
            r = ones(40,1)*obj.barRad;
            [xx,yy,zz] = cylinder(r);
            obj.barTForm = gobjects(obj.bb,1);
            for i = 1:obj.bb
                barHandle = surf(ax,xx,yy,zz,'LineStyle', 'none');
                obj.barTForm(i) = hgtransform('Parent',ax);
                set(barHandle,'Parent',obj.barTForm(i))
                %updateMember(obj.nodePoints(obj.barNodes(1,i),:),obj.nodePoints(obj.barNodes(2,i),:),...
                %   obj.barTForm(i));
            end
            
            r = ones(40,1)*obj.stringRad;
            [xx,yy,zz] = cylinder(r);
            obj.stringTForm = gobjects(obj.ss,1);
            for i = 1:obj.ss
                stringHandle = surf(ax,xx,yy,zz+100,'LineStyle', 'none');
                obj.stringTForm(i) = hgtransform('Parent',ax);
                stringTFormOffset = hgtransform('Parent',obj.stringTForm(i));
                set(stringHandle,'Parent',stringTFormOffset)
                set(stringTFormOffset,'matrix',makehgtform('translate',[0 0 -100]));
                %updateMember(obj.nodePoints(obj.stringNodes(1,i),:),obj.nodePoints(obj.stringNodes(2,i),:),...
                %   obj.stringTForm(i));
            end
        end
        function updatePlot(obj)
            H = eye(4);
            for i = 1:obj.n
                H(1:3,4) = obj.nodePoints(i,:)';
                obj.sphereTForm(i).Matrix = H;
            end
           % H = zeros(obj.bb+obj.ss,4,4);
           if(isempty(obj.barNodes))
               nodeXYZ1 = obj.nodePoints( obj.stringNodes(1,:),:);
            nodeXYZ2 = obj.nodePoints( obj.stringNodes(2,:),:);
           else
            nodeXYZ1 = obj.nodePoints([obj.barNodes(1,:), obj.stringNodes(1,:)],:);
            nodeXYZ2 = obj.nodePoints([obj.barNodes(2,:), obj.stringNodes(2,:)],:);
           end          
            tForms = [obj.barTForm; obj.stringTForm];
            updateMember(nodeXYZ1,nodeXYZ2,tForms);
        end
        
    end
end

function updateMember(p1,p2,tform)
length = sum((p2-p1).^2,2).^0.5;
n = size(p1,1);
v1 = repmat([0 0 1],n,1);
v2 = (p2-p1)./repmat(length,1,3);
d = sum(v1 .* v2,2);
isVert= 1:n;
isVert = isVert(d < -0.999999999);
isOtherVert= 1:n;
isOtherVert = isOtherVert(d > 0.9999999999);
s = ((1+d)*2).^0.5;
axis = cross(v1,v2);
if(~isempty(isOtherVert))
    axis(isOtherVert,:) = repmat([0.00001 0 0],size(isOtherVert,2),1);
    s(isOtherVert) = ones(size(isOtherVert,2),1);
end
if(~isempty(isVert))
    axis(isVert,:) = repmat([1 0 0],size(isVert,2),1);
    s(isVert) = 0.0001*ones(size(isVert,2),1);
end
q = [0.5*s, axis./repmat(s,1,3)];
H1  = repmat(eye(4),1,1,n);
H1(1:3,1:3,:) = permute(quat2dcm(q),[2 1 3]);
H1(1:3,3,:) = H1(1:3,3,:).*repmat(reshape(length,1,1,[]),3,1,1); 
H1(1:3,4,:) = reshape(p1',3,1,[]);
for i = 1:n
    tform(i).Matrix = H1(:,:,i);
end
end
