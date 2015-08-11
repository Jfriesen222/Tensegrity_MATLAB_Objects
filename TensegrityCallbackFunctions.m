classdef TensegrityCallbackFunctions <handle
    properties
        vals % vector of values which can be updated 
    end
    
    methods
        function obj = TensegrityCallbackFunctions(vals)
            obj.vals = vals; %initialize values
        end
        
        function updateVal(obj,val,i)
            obj.vals(i) = val; %update vals at the ith index
        end
        
        function timerUpdate(obj,funcHandle)
            %this is a function wrapper for an update function
            vec = obj.vals;
            funcHandle(vec);
        end
    end
    methods(Static)
        function makeSlider(f, x, y,func,lims,title)
            %Fuction to make a quick slider
            %f - a handle to the figure the slider will be created in
            %x - the x position in pixels of the slider
            %y - the y position in pixels of the slider
            %func - the callback function the slider will execute when
            %       moved
            %lims - lower and upper limits for the slider in a length 2
            %       vector
            %title - a string which contains the title for the slider
            
            b2 = uicontrol('Parent',f,'Style','slider','Position',[51,54,419,23]+[x y 0 0],...
                'value',mean(lims), 'min',lims(1), 'max',lims(2));
            bgcolor = f.Color;
            uicontrol('Parent',f,'Style','text','Position',[45,25,23,23]+[x y 0 0],...
                'String',num2str(lims(1),3),'BackgroundColor',bgcolor);
            uicontrol('Parent',f,'Style','text','Position',[450,25,23,23]+[x y 0 0],...
                'String',num2str(lims(2),3),'BackgroundColor',bgcolor);
            uicontrol('Parent',f,'Style','text','Position',[210,25,100,23]+[x y 0 0],...
                'String',title,'BackgroundColor',bgcolor);
            
            %only way I found to get the slider to update smoothly when
            %dragged
            addlistener(b2,'ContinuousValueChange', func);
        end
    end
    
end

