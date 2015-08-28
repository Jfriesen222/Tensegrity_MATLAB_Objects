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
        memberTForms              %transform object for bar cylinders
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
            [xx,yy,zz] = sphere(9); xx = xx*obj.barRad; yy = yy*obj.barRad; zz = zz*obj.barRad;
            obj.sphereTForm = gobjects(obj.n,1);
            hold(ax,'on');
            for i = 1:obj.n
                sphereHandle = surf(ax,xx,yy,zz,'LineStyle', 'none');
                obj.sphereTForm(i) = hgtransform('Parent',ax);
                set(sphereHandle,'Parent',obj.sphereTForm(i))
                set(obj.sphereTForm(i),'matrix',makehgtform('translate',obj.nodePoints(i,:)));
            end
            cylRes = 20;
            r = ones(2,1)*obj.barRad;
            [xx,yy,zz] = cylinder(r,cylRes);
            barTForm = gobjects(obj.bb,1);
            for i = 1:obj.bb
                barHandle = surf(ax,xx,yy,zz,'LineStyle', 'none');
                barTForm(i) = hgtransform('Parent',ax);
                set(barHandle,'Parent',barTForm(i))
            end
            
            r = ones(2,1)*obj.stringRad;
            [xx,yy,zz] = cylinder(r,cylRes);
            stringTForm = gobjects(obj.ss,1);
            for i = 1:obj.ss
                stringHandle = surf(ax,xx,yy,zz+100,'LineStyle', 'none');
                stringTForm(i) = hgtransform('Parent',ax);
                stringTFormOffset = hgtransform('Parent',stringTForm(i));
                set(stringHandle,'Parent',stringTFormOffset)
                set(stringTFormOffset,'matrix',makehgtform('translate',[0 0 -100]));
            end
            obj.memberTForms = [barTForm; stringTForm];
        end
        function updatePlot(obj)
            nn = obj.n;
            HH = cell(nn,1);
            for i =1:nn
                HH{i} = [eye(3), obj.nodePoints(i,:)'; 0 0 0 1] ;
            end
           [obj.sphereTForm(:).Matrix] = HH{:};
            if(isempty(obj.barNodes))
                nodeXYZ1 = obj.nodePoints( obj.stringNodes(1,:),:);
                nodeXYZ2 = obj.nodePoints( obj.stringNodes(2,:),:);
            else
                nodeXYZ1 = obj.nodePoints([obj.barNodes(1,:), obj.stringNodes(1,:)],:);
                nodeXYZ2 = obj.nodePoints([obj.barNodes(2,:), obj.stringNodes(2,:)],:);
            end
            updateMember(nodeXYZ1,nodeXYZ2,obj.memberTForms);
        end
        
    end
end

function updateMember(p1,p2,tform)
n = size(p1,1);
length = sum((p2-p1).^2,2).^0.5;
lengthInv = ones(n,1)./length;
vec = (p2-p1).*lengthInv(:,[1 1 1]);
d =  vec(:,3);
isVert= 1:n;
isVert = isVert(d < -0.999999999);
isOtherVert= 1:n;
isOtherVert = isOtherVert(d > 0.9999999999);
s = ((1+d)*2).^0.5;
axis = vec(:,[2 1]);
axis(:,1) = -axis(:,1);
if(~isempty(isOtherVert))
    axis(isOtherVert,:) = repmat(sin(pi-0.001)*[1 0],size(isOtherVert,2),1);
    s(isOtherVert) = ones(size(isOtherVert,2),1);
end
if(~isempty(isVert))
    axis(isVert,:) = repmat([1 0],size(isVert,2),1);
    s(isVert) = 0.0000001*ones(size(isVert,2),1);
end
qin = [0.5*s, axis./s(:,[1 1])];
qinmag = sqrt(sum(qin.^2,2));
qin = qin./qinmag(:,[1 1 1]);
H1  = repmat(eye(4),1,1,n);
H1(1,1,:) = qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2;
H1(2,1,:) = 2.*(qin(:,2).*qin(:,3));
H1(3,1,:) = 2.*( - qin(:,1).*qin(:,3));
H1(1,2,:) = 2.*(qin(:,2).*qin(:,3)) ;
H1(2,2,:) = qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2;
H1(3,2,:) = 2.*( + qin(:,1).*qin(:,2));
H1(1,3,:) = 2*length.*( + qin(:,1).*qin(:,3));
H1(2,3,:) = 2*length.*( - qin(:,1).*qin(:,2));
H1(3,3,:) = length.*( qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2);
H1(1:3,4,:) = reshape(p1',3,1,[]);
HH = cell(n,1);
for i =1:n
    HH{i} = H1(:,:,i);
end
[tform(:).Matrix] = HH{:};
end
