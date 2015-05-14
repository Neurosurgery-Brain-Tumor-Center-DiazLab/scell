classdef ScatterSelect < handle & MObject
  %SCATTERSELECT Interactively select points from a scatterplot
  %   Detailed explanation goes here
  
properties
  data % Nx2 matrix of data points, each row = [x y]
  title = 'Unnamed' % figHure title
  closeFcn % function handle is called when the figHure is closed
  symbol = 'O' % current plot symbol  
  highColor = 'red' % highligh color
end

properties (GetAccess = public, SetAccess = private)
  figH     % handle to figure
  figHAxH  % handle to figure axes
  scatterH % handle to scatter object
  highH % handle to higlight object
  highInd % index to highlight value
  selecting = false
  isInWindow = false % whether mouse pointer is inside the window
  isInTimer % timer for checking if mouse pointer is inside the window
  lastSettings
end

methods
  function self = ScatterSelect()
    self@MObject();
    self.isInTimer = timer('TimerFcn', @self.checkIsIn, 'Period', 0.1,...
      'ExecutionMode', 'fixedSpacing');
    self.lastSettings = self.defaultSettings;
  end

  function set.title(self, value)
    self.title = value;
    self.updatePlot();
  end

  function set.data(self, value)
    self.data = value;
    self.updatePlot();
  end

  function val = defaultSettings(self)
    val = [86 344 672  504];
  end
  
  function val = getSettings(self)
    if ishandle(self.figH)
      val = get(self.figH, 'Position');
    else
      val = self.defaultSettings();
    end
    self.lastSettings = val;
  end
  
  function changeToSettings(self, settings)
    self.lastSettings = settings;
    self.updatePlotToSettings(settings);
  end
  
  function show(self)
    if ishandle(self.figH)
      set(self.figH, 'Visible', 'on');
      stop(self.isInTimer);
      start(self.isInTimer);
    else
      f = figure();
%       xlim([0 2]);
%       ylim([0 3]); % debugging
      set(f, 'Units', 'pixels', 'CloseRequestFcn', ...
        @self.windowAboutToClose, 'WindowButtonMotionFcn', ...
        @self.mouseMoved);        
      self.figH = f;
      self.figHAxH = gca;
      % Keeping the units normalized allows axes to auto resize,
      % strange side effect of setting a property.
      set(self.figHAxH, 'Unit', 'normalized');
      start(self.isInTimer);
      self.updatePlotToSettings(self.lastSettings);
    end            
    self.updatePlot();
  end

  function selectBegin(self)
    self.selecting = true;
  end

  function selectEnd(self)
    self.selecting = false;
  end
  
  function closeFigure(self)
  % Closes figure without invoking the closing callback
    if ishandle(self.figH)
      set(self.figH, 'CloseRequestFcn', '');
      stop(self.isInTimer);
      delete(self.figH);
    end
  end    
  
end

methods (Access = private)
  function windowAboutToClose(self, varargin)
    stop(self.isInTimer);    
    if isa(self.closeFcn, 'function_handle')
      self.closeFcn();
    else 
      delete(self.figH); % debugging
    end
  end  
  
  function mouseMoved(self, varargin)
    if ~isempty(self.data)
      % find
      ind = self.closestDataToPointer();
      delete(self.highH);
      hold on;
      self.highH = plot(self.figHAxH, self.data(ind,1), self.data(ind,2), ...
        'Marker', self.symbol, 'Color', self.highColor);
      hold off;
    end
  end
  
  function ind = closestDataToPointer(self)
  % find index to the data point closest to the current mouse pointer
  % returns empty if not found
    ind = [];
    if ~isempty(self.data)
      xy = self.pointerLocInCoord();
      vecdiff = self.data - repmat(xy, [length(self.data) 1]);
      diff = zeros(1,length(vecdiff));
      for i=1:length(diff)
        diff(i) = norm(vecdiff(i,:));
      end
      [~, ind] = min(diff);
    end
  end
  
  function checkIsIn(self, varargin)
    fp = get(self.figH, 'Position');
    pointer = get(groot, 'PointerLocation');
    isIn = false;
    if self.pointInRect(fp, pointer(1), pointer(2))
      isIn = true;
    end
    if ~isIn
      delete(self.highH);
%       set(self.figH, 'Name', 'Mouse OUT'); %debugging
    else
%       set(self.figH, 'Name', 'Mouse IN'); %debugging
    end
    self.isInWindow = isIn;    
%     [x,y]=self.pointerLocInCoord(); % debugging
%     fprintf('(%f, %f)\n', x, y);    
  end
  
  function updatePlot(self)
    if ishandle(self.figH)
      set(self.figH, 'Name', self.title);
      if isempty(self.data)
        clf(self.figH);
      else
        self.scatterH = scatter(self.figHAxH, self.data(:,1), self.data(:,2));
      end
      drawnow;
    end    
  end  
  
  function updatePlotToSettings(self, settings)
    if ishandle(self.figH)
      set(self.figH, 'Position', settings);
    end
  end

  function xy = pointerLocInCoord(self)
    pl = get(groot, 'PointerLocation');
    fp = get(self.figH, 'Position');
    ap = get_in_units(self.figHAxH, 'Position', 'Pixels');
    xl = xlim(self.figHAxH);
    yl = ylim(self.figHAxH);
    pixToCoordXScale = (xl(2) - xl(1)) / ap(3);
    pixToCoordYScale = (yl(2) - yl(1)) / ap(4);    
    x = (pl(1) - (fp(1) + ap(1))) * pixToCoordXScale;
    y = (pl(2) - (fp(2) + ap(2))) * pixToCoordYScale;
    xy = [x y];
  end
  
  function tf = pointInRect(self, rect, x,y)
  % rect = [x y width height]
  % x,y test point
    tf = false;
    if (rect(1) <= x) && (x <= rect(1) + rect(3)) && ...
      (rect(2) <= y) && (y <= rect(2) + rect(4))
      tf = true;
    end
  end
end
  
end

