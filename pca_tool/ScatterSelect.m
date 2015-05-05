classdef ScatterSelect < handle
  %SCATTERSELECT Interactively select points from scatterplot
  %   Detailed explanation goes here
  
properties
  data % Nx2 matrix of data points, each row = [x y]
  title = 'Unnamed' % figure title
  closeFcn % function handle is called when the figure is closed
  symbol = 'O' % current plot symbol  
  highColor = 'red' % highligh color
end

properties (GetAccess = public, SetAccess = private)
  fig
  figAx  
  scatterH % handle to scatter object
  highH % handle to higlight object
  highInd % index to highlight value
  selecting = false
  isInWindow = false % whether mouse pointer is inside the window
  isInTimer % timer for checking if mouse pointer is inside the window
  pixToCoordXScale % scaling from pixels to figure coordinates, x
  pixToCoordYScale % scaling from pixels to figure coordinates, y
end

methods
  function self = ScatterSelect()
    self.isInTimer = timer('TimerFcn', @self.checkIsIn, 'Period', 0.1,...
      'ExecutionMode', 'fixedSpacing');
  end

  function set.title(self, value)
    self.title = value;
    self.updatePlot();
  end

  function set.data(self, value)
    self.data = value;
    self.updatePlot();
  end

  function show(self)
    if ishandle(self.fig)
      set(self.fig, 'Visible', 'on');
      stop(self.isInTimer);
      start(self.isInTimer);
    else
      f = figure();
%       xlim([0 2]);
%       ylim([0 3]); % debugging
      set(f, 'Units', 'pixels', 'CloseRequestFcn', ...
        @self.windowAboutToClose, 'WindowButtonMotionFcn', ...
        @self.mouseMoved);        
      self.fig = f;
      self.figAx = gca;
      set(self.figAx, 'Unit', 'pixels');
      start(self.isInTimer);
    end            
    self.updatePlot();
  end

  function selectBegin(self)
    self.selecting = true;
  end

  function selectEnd(self)
    self.selecting = false;
  end
  
end

methods (Access = private)
  function windowAboutToClose(self, varargin)
    stop(self.isInTimer);
    delete(self.fig); % debugging
    if isa(self.closeFcn, 'function_handle')
      self.closeFcn();
    end
  end
  
  function updatePixToCoordScaling(self)    
    xl = xlim(self.figAx);
    yl = ylim(self.figAx);
    ap = get(self.figAx, 'Position');
    self.pixToCoordXScale = (xl(2) - xl(1)) / ap(3);
    self.pixToCoordYScale = (yl(2) - yl(1)) / ap(4);
  end
  
  function mouseMoved(self, varargin)
    % find
    ind = self.closestDataToPointer();
    delete(self.highH);
    hold on;
    self.highH = plot(self.figAx, self.data(ind,1), self.data(ind,2), ...
      'Marker', self.symbol, 'Color', self.highColor);
    hold off;
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
    fp = get(self.fig, 'Position');
    pointer = get(groot, 'PointerLocation');
    isIn = false;
    if self.pointInRect(fp, pointer(1), pointer(2))
      isIn = true;
    end
    if ~isIn
      delete(self.highH);
      set(self.fig, 'Name', 'Mouse OUT'); %debugging
    else
      set(self.fig, 'Name', 'Mouse IN'); %debugging
    end
    self.isInWindow = isIn;    
%     [x,y]=self.pointerLocInCoord(); % debugging
%     fprintf('(%f, %f)\n', x, y);    
  end
  
  function updatePlot(self)
    if ishandle(self.fig)
      set(self.fig, 'Name', self.title);
      if isempty(self.data)
        clf(self.fig);
      else
        self.scatterH = scatter(self.figAx, self.data(:,1), self.data(:,2));
      end
      drawnow;
      self.updatePixToCoordScaling();
    end    
  end  

  function xy = pointerLocInCoord(self)
    pl = get(groot, 'PointerLocation');
    fp = get(self.fig, 'Position');
    ap = get(self.figAx, 'Position');
    x = (pl(1) - (fp(1) + ap(1))) * self.pixToCoordXScale;
    y = (pl(2) - (fp(2) + ap(2))) * self.pixToCoordYScale;
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

