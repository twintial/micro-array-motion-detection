classdef (StrictDefaults)MyDoaDisplay < matlab.System
    %DOADisplay Helper object for Direction-Of-Arrival display
    %   H = DOADisplay returns a DOA display.
    %
    %   H = DOADisplay('Name', Value, ...) returns a DOA
    %   display with each specified property name set to the specified
    %   value. You can specify additional name-value pair arguments in any
    %   order as (Name1,Value1,...,NameN,ValueN).
    %
    %   Step method syntax:
    %
    %   step(H, X) plots the subsequent directions in vector X one after
    %   the other.
    %
    %   System objects may be called directly like a function instead of
    %   using the step method. For example, y = step(obj, x) and y = obj(x)
    %   are equivalent.
    %
    %   DOADisplay methods:
    %
    %   step               - See above description for use of this method
    %   release            - Allow property value and input characteristics
    %                        changes
    %   delete             - Closes display window
    %
    %   DOADisplay properties:
    %
    %   TipRadius           - Length of the arrow, from origin to tip
    %   SidesRadius         - Distance on the ends of the arrow sides from
    %                         the origin
    %   SidesAngleInRadians - Angular span of the arrow sides from the
    %                         arrow tip
    %
    %
    %   This class DOADisplay is only in support of
    %   AudioArrayDirectionOfArrivalEstimationExample. It may change in a
    %   future release.
    
    % Copyright 2013-2018 The MathWorks, Inc.
    
    properties(Access = private)
        ArrowLineHandles
    end
    
    properties(Access = public)
        TipRadius = 1
        SidesRadius = 0.9
        SidesAngleInRadians = 0.05
    end
    
    methods(Access = protected)
        function setupImpl(obj, ~)
            initFigure(obj);
        end
        function stepImpl(obj, theta_rad)
            if ~ishandle(obj.ArrowLineHandles)
                initFigure(obj);
            end
            for k = 1:numel(theta_rad)
                if ishandle(obj.ArrowLineHandles)
                    [xdata, ydata] = getArrowXYData(obj, theta_rad(k));
                    set(obj.ArrowLineHandles, 'XData', xdata, 'YData', ydata)
                    drawnow limitrate;
                end
            end

        end
    end
    
    methods
        function obj = DOADisplay(varargin)
            %Constructor
            
            % Support name-value pair arguments
            setProperties(obj, nargin, varargin{:});
        end
        function delete(obj)
            delete@matlab.System(obj);
            if ishandle(obj.ArrowLineHandles)
                fig = get(get(obj.ArrowLineHandles,'Parent'),'Parent');
                close(fig)
            end
            
        end
    end
    
    methods(Access = private)
        function initFigure(obj)
            figure
            obj.ArrowLineHandles = compass(1, 0);
            [x0, y0] = getArrowXYData(obj, 0);
            set(obj.ArrowLineHandles, ... 
                'XData', x0, 'YData', y0, ...
                'LineWidth', 3, ...
                'LineStyle','-',...
                'Color','r')
            
            % Only get half polar plot
            % set(get(obj.ArrowLineHandles,'Parent'), 'Xlim', [0,1])
            txt = findall(get(obj.ArrowLineHandles,'Parent'), 'Type', 'Text');
            unwantedText = cellstr(num2str((120:30:240)'));
            toRemove = false(size(txt));
            for kr = 1:length(unwantedText)
                toRemove = toRemove|strcmp(get(txt, 'String'), unwantedText{kr});
            end
            set(txt(toRemove),'Visible','off')
        end
        function [xdata, ydata] = getArrowXYData(obj, theta_rad)

            xtip = obj.TipRadius * cos(theta_rad);
            ytip = obj.TipRadius * sin(theta_rad);

            xleftside = obj.SidesRadius * cos(theta_rad + obj.SidesAngleInRadians);
            yleftside = obj.SidesRadius * sin(theta_rad + obj.SidesAngleInRadians);

            xrightside = obj.SidesRadius * cos(theta_rad - obj.SidesAngleInRadians);
            yrightside = obj.SidesRadius * sin(theta_rad - obj.SidesAngleInRadians);

            xdata = [0 xtip xleftside xtip xrightside];
            ydata = [0 ytip yleftside ytip yrightside];


        end
    end
end