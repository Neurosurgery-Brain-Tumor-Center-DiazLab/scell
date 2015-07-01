classdef ScatterSelect < MObject
  %SCATTERSELECT Interactively select points from a scatterplot
  %   Detailed explanation goes here
  
properties
  title = 'Unnamed' % figure title
  closeFcn % function handle is called when the figure is closed
  marker = 'O' % current plot marker
  color = [0 0.4470 0.7410] % default Matlab color
  highColor = 'red' % highligh color
  highMarker = '*'  
  selectingEnabled = true
end

properties (GetAccess = public, SetAccess = private)
%   data % Nx2 matrix of data points, each row = [x y]
%   cluster % empty or Nx1 vector of cluster tags (integers >= 1), 
%           % each corresponds to data row
%   selData % selected data Mx2
%   selIndices % incides to self.data corresponding to the selected data  
%   currentIndex % index to data closest to the mouse button now, empty if
%                % mouse is not inside the figure
  pm % ScatterSelectPM presentation model
  figH     % handle to figure
  figAxH  % handle to figure axes
  scatterH % handle to scatter object
  highH % handle to higlight object
  highInd % index to highlight value
  selH % handle to selected points
  isInWindow = false % whether mouse pointer is inside the window
  isInTimer % timer for checking if mouse pointer is inside the window
  lastSettings
%   clusterCount = 1 % number of clusters
  colorMap = brewermap(8,'Dark2');
  plotBorder = 0.05; % fraction of free space (border) around xy plot
end


methods
  function self = ScatterSelect()
    self@MObject();
    s = ScatterSelectPM();
    s.connectMe('selection_changed', @self.updatePlotSelection);
    s.connectMe('data_changed', @self.updatePlot);
    self.pm = s;
    self.lastSettings = self.defaultSettings;    
  end

  function set.title(self, value)
    self.title = value;
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
      f = figure('HandleVisibility', 'off');
      self.figAxH = axes('Parent', f);
%       xlim([0 2]);
%       ylim([0 3]); % debugging
      set(f, 'Units', 'pixels', 'CloseRequestFcn', ...
        @self.windowAboutToClose, 'WindowButtonMotionFcn', ...
        @self.mouseMoved, 'WindowButtonDownFcn', @self.mouseClicked, ...
        'KeyPressFcn', @self.keyPressed);        
      self.figH = f;      
      hold(self.figAxH, 'on');
      self.selH = plot(0, 0, 'Visible', 'off', 'Parent', self.figAxH,...
        'Color', self.highColor, 'Marker', self.highMarker, ...
        'LineStyle', 'none');      
      self.highH = plot(0, 0, 'Visible', 'off', 'Parent', self.figAxH,...
        'Color', self.highColor, 'Marker', self.highMarker);
      self.scatterH = scatter(0, 0, 'Visible', 'off', 'Parent', ...
        self.figAxH, 'Marker', self.marker);
      hold(self.figAxH, 'off');
      % Keeping the units normalized allows axes to auto resize,
      % strange side effect of setting a property.
      set(self.figAxH, 'Unit', 'normalized');
      self.isInTimer = timer('TimerFcn', @self.checkIsIn, 'Period', 0.1,...
      'ExecutionMode', 'fixedSpacing', 'Name', 'is_in_timer');      
      start(self.isInTimer);
      self.updatePlotToSettings(self.lastSettings);
    end            
    self.updatePlot();
  end
  
  function closeFigure(self)
  % Closes figure without invoking the closing callback
    if ishandle(self.figH)
      set(self.figH, 'CloseRequestFcn', '');
      if isvalid(self.isInTimer)
        stop(self.isInTimer);
        delete(self.isInTimer);
        clear self.isInTimer;
      end
      delete(self.figH);
    end
  end    
    
end

%*** Private 
methods (Access = private)
  function windowAboutToClose(self, varargin)
    if isvalid(self.isInTimer)
      stop(self.isInTimer);  
      delete(self.isInTimer);
      clear self.isInTimer;
    end
    if isa(self.closeFcn, 'function_handle')      
      self.closeFcn();
    else 
      if ishandle(self.figH)
        delete(self.figH); % debugging
      end
    end
  end  
  
  function mouseMoved(self, varargin)
    if ~isempty(self.pm.data) && self.isInWindow
      % find
      ind = self.closestDataToPointer();
      set(self.highH, 'XData', self.pm.data(ind,1), 'YData', ...
        self.pm.data(ind,2), 'Visible', 'on', 'Color', self.highColor);
%       self.emit('highlight', ind);
      self.pm.updateCurrentIndex(ind);
    end
  end
  
  function mouseClicked(self, varargin)
    if self.selectingEnabled
      self.pm.selectCurrentIndex();
    end
  end
  
  function keyPressed(self, varargin)
    ch = get(self.figH, 'CurrentCharacter');
    isEsc = false;
    isClear = false;
    % capture buttons
    if ~isempty(ch)
      if uint16(ch) == 27
        isEsc = true;
      elseif lower(ch) == 'c'
        isClear = true;
      end
    end
    % actions
    if self.selectingEnabled
      if isEsc
        self.pm.deselectLastIndex();
      elseif isClear
        self.pm.deselectAll();
      end
    end
  end
  
  function ind = closestDataToPointer(self)
  % find index to the data point closest to the current mouse pointer
  % returns empty if not found
    ind = [];
    if ~isempty(self.pm.data)
      xy = self.pointerLocInCoord();
      vecdiff = self.pm.data - repmat(xy, [length(self.pm.data) 1]);
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
      set(self.highH, 'Visible', 'off');
      self.pm.updateCurrentIndex([]);
    end    
    if self.isInWindow ~= isIn      
      self.emit('is_in', isIn);
    end
    self.isInWindow = isIn;
  end
  
  function updatePlot(self)
    if self.pm.clusterCount > size(self.colorMap,1)
      warning('Specified more clusters than there are different colors');
    end
    if ishandle(self.figH)
      set(self.figH, 'Name', self.title); 
      if ~isempty(self.pm.data)
%         set(self.figAxH, 'XLimMode', 'auto', 'YLimMode', 'auto');   
        cmap = self.computeClusterColors();
        set(self.scatterH, 'XData', self.pm.data(:,1), 'YData', ...
          self.pm.data(:,2), 'Visible', 'on', 'CData', cmap);
        minX = min(self.pm.data(:,1));
        maxX = max(self.pm.data(:,1));
        dx = maxX(1) - minX(1);
        minY = min(self.pm.data(:,2));
        maxY = max(self.pm.data(:,2));
        dy = maxY(1) - minY(1);
        offX = self.plotBorder * dx;
        offY = self.plotBorder * dy;
        set(self.figAxH, 'XLimMode', 'manual', 'YLimMode', 'manual',...
          'XLim', [minX(1)-offX maxX(1)+offX], 'YLim', ...
          [minY(1)-offY maxY(1)+offY]);
      end
      self.updatePlotSelection();
      drawnow;
    end    
  end  
  
  function updatePlotSelection(self)
    if ~isempty(self.pm.selData)
      set(self.selH, 'XData', self.pm.selData(:,1), 'YData', ...
        self.pm.selData(:,2), 'Visible', 'on');      
    else
      set(self.selH, 'Visible', 'off');
    end
  end
  
  function cmap = computeClusterColors(self)
    if self.pm.clusterCount > 1      
      cmap = zeros(size(self.pm.data,1), 3);
      colorCount = size(self.colorMap,1);
      for i = 1:size(self.pm.data,1)
        ind = mod(self.pm.cluster(i), colorCount) + 1;
        cmap(i,:) = self.colorMap(ind,:);
      end      
    else
      cmap = self.color;    
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
    ap = get_in_units(self.figAxH, 'Position', 'Pixels');
    xl = xlim(self.figAxH);
    yl = ylim(self.figAxH);
    pixToCoordXScale = (xl(2) - xl(1)) / ap(3);
    pixToCoordYScale = (yl(2) - yl(1)) / ap(4);    
    x = (pl(1) - (fp(1) + ap(1))) * pixToCoordXScale;
    y = (pl(2) - (fp(2) + ap(2))) * pixToCoordYScale;
    xy = [x y] + [xl(1) yl(1)];
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

