classdef ScatterSelect < handle & MObject
  %SCATTERSELECT Interactively select points from a scatterplot
  %   Detailed explanation goes here
  
properties
  title = 'Unnamed' % figHure title
  closeFcn % function handle is called when the figHure is closed
  marker = 'O' % current plot marker
  color = [0 0.4470 0.7410] % default Matlab color
  highColor = 'red' % highligh color
  highMarker = '*'  
  selectingEnabled = true
end

properties (GetAccess = public, SetAccess = private)
  data % Nx2 matrix of data points, each row = [x y]
  cluster % empty or Nx1 vector of cluster tags (integers >= 1), 
          % each corresponds to data row
  selData % selected data Mx2
  selIndices % incides to self.data corresponding to the selected data  
  currentIndex % index to data closest to the mouse button now, empty if
               % mouse is not inside the figure
  figH     % handle to figure
  figAxH  % handle to figure axes
  scatterH % handle to scatter object
  highH % handle to higlight object
  highInd % index to highlight value
  selH % handle to selected points
  isInWindow = false % whether mouse pointer is inside the window
  isInTimer % timer for checking if mouse pointer is inside the window
  lastSettings
  clusterCount = 1 % number of clusters
  colorMap = brewermap(8,'Dark2');
  uiState % current ui state
end


methods
  function self = ScatterSelect()
    self@MObject();
    self.lastSettings = self.defaultSettings;    
  end

  function set.title(self, value)
    self.title = value;
    self.updatePlot();
  end

%   function set.data(self, value)
%     self.data = value;
%     self.updatePlot();
%   end

  function updateData(self, data, cluster)
    self.data = data;
    self.clusterCount = 1;
    self.cluster = cluster;
    if ~isempty(cluster)
      tags = [];
      for i=1:length(cluster)
        if ~ismember(cluster(i), tags)
          tags(end+1) = cluster(i); %#ok<AGROW>
        end
      end
      self.clusterCount = length(tags);
    end
    if self.clusterCount > size(self.colorMap,1)
      warning('Specified more clusters than there are different colors');
    end
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
      'ExecutionMode', 'fixedSpacing');      
      start(self.isInTimer);
      self.updatePlotToSettings(self.lastSettings);
    end            
    self.updatePlot();
  end

  function clearSelection(self)
    self.selData = [];
    self.selIndices = [];
    set(self.selH, 'Visible', 'off', 'XData', [], 'YData', []);
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
    if ~isempty(self.data) && self.isInWindow
      % find
      ind = self.closestDataToPointer();
      set(self.highH, 'XData', self.data(ind,1), 'YData', ...
        self.data(ind,2), 'Visible', 'on', 'Color', self.highColor);
      self.emit('highlight', ind);
      self.currentIndex = ind;
    end
  end
  
  function mouseClicked(self, varargin)
    if self.selectingEnabled
      if ~ismember(self.currentIndex, self.selIndices)
        self.selIndices(end+1) = self.currentIndex;
        self.selData = self.data(self.selIndices, :);
        self.updatePlotSelection();
        self.emit('selection', self.selData, self.selIndices);
      end
    end
  end
  
  function keyPressed(self, varargin)
    ch = get(self.figH, 'CurrentCharacter');
    if ~isempty(ch) && (uint16(ch) == 27)
      isEsc = true;
    else
      isEsc = false;
    end
    % escape pressed ?
    if self.selectingEnabled && isEsc
      self.selIndices = self.selIndices(1:end-1);
      self.selData = self.data(self.selIndices, :);
      self.updatePlotSelection();
      self.emit('selection', self.selData, self.selIndices);
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
      set(self.highH, 'Visible', 'off');
      self.currentIndex = [];
    end    
    if self.isInWindow ~= isIn
      self.emit('is_in', isIn);
    end
    self.isInWindow = isIn;
  end
  
  function updatePlot(self)
    if ishandle(self.figH)
      set(self.figH, 'Name', self.title); 
      if ~isempty(self.data)
        set(self.figAxH, 'XLimMode', 'auto', 'YLimMode', 'auto');   
        cmap = self.computeClusterColors();
        set(self.scatterH, 'XData', self.data(:,1), 'YData', ...
          self.data(:,2), 'Visible', 'on', 'CData', cmap);
%         self.scatterH = scatter(self.figAxH, self.data(:,1), ...
%           self.data(:,2));
        set(self.figAxH, 'XLimMode', 'manual', 'YLimMode', 'manual');
      end
      self.updatePlotSelection();
    end    
  end  
  
  function updatePlotSelection(self)
    if ~isempty(self.selData)
      set(self.selH, 'XData', self.selData(:,1), 'YData', ...
        self.selData(:,2), 'Visible', 'on');      
    else
      set(self.selH, 'Visible', 'off');
    end
  end
  
  function cmap = computeClusterColors(self)
    if self.clusterCount > 1      
      cmap = zeros(size(self.data,1), 3);
      colorCount = size(self.colorMap,1);
      for i = 1:size(self.data,1)
        ind = mod(self.cluster(i), colorCount) + 1;
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

