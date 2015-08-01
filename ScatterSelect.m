classdef ScatterSelect < MObject
  %SCATTERSELECT Interactively select points from a scatterplot
  %   Detailed explanation goes here
  
properties  
  closeFcn % function handle is called when the figure is closed
  marker = 'O' % current plot marker
  color = [0 0.4470 0.7410] % default Matlab color
  highColor = 'red' % highligh color
  highMarker = '*'    
  selectingEnabled = true
end

properties (GetAccess = public, SetAccess = private)
  pm % ScatterSelectPM presentation model
  figH     % handle to figure
  figAxH  % handle to figure axes
  scatterH % handle to scatter object
  highH % handle to higlight object
  highInd % index to highlight value
  selH % handle to selected points
  lassoH % handle to lasso select polygon
  isInWindow = false % whether mouse pointer is inside the window
  isInTimer % timer for checking if mouse pointer is inside the window
  lastSettings
  cmap
  clusterCount = 1 % number of clusters
  colorMap = brewermap(8,'Dark2');
  plotBorder = 0.05; % fraction of free space (border) around xy plot
end


methods
  function self = ScatterSelect()
    self@MObject();
    s = ScatterSelectPM();
    s.connectMe('selection_changed', @self.updatePlotSelection);
    s.connectMe('data_changed', @self.updatePlot);
    s.connectMe('changed', @self.refresh);
    self.pm = s;
    self.lastSettings = self.defaultSettings;    
  end

%   function set.title(self, value)
%     self.title = value;
%     self.updatePlot();
%   end
  
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
      %*** Start creating UI related graphics elements
      hold(self.figAxH, 'on');
      % selected points
      self.selH = plot(0, 0, 'Visible', 'off', 'Parent', self.figAxH,...
        'Color', self.highColor, 'Marker', self.highMarker, ...
        'LineStyle', 'none');      
      % currently highlighted point
      self.highH = plot(0, 0, 'Visible', 'off', 'Parent', self.figAxH,...
        'Color', self.highColor, 'Marker', self.highMarker);
      % lasso select polygon
      self.lassoH = line('Visible', 'off', 'Parent', self.figAxH,...
        'Color', self.pm.lassoColor);
      % the actual data points
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
      if self.pm.selectionMode == SelectionMode.Lasso
        self.refreshLasso();
      end
    end
  end
  
  function mouseClicked(self, varargin)
    if self.selectingEnabled
      switch self.pm.selectionMode
        case SelectionMode.Pointwise
          self.pm.selectCurrentIndex();
        case SelectionMode.Lasso
          self.pm.addLassoPoint(self.pointerLocInCoord());
        otherwise
          error('Bug found!');
      end
      
    end
  end
  
  function keyPressed(self, varargin)
    ch = get(self.figH, 'CurrentCharacter');
    isEsc = false;
    isC = false;
    isQ = false;
    is1 = false;
    is2 = false;
    isEnter = false;
    % capture buttons
    if ~isempty(ch)
      if uint16(ch) == 27
        isEsc = true;
      elseif uint16(ch) == 13
        isEnter = true;
      elseif lower(ch) == 'c'
        isC = true;
      elseif lower(ch) == 'q'
        isQ = true;
      elseif lower(ch) == '1'
        is1 = true;
      elseif lower(ch) == '2'
        is2 = true;
      end
    end
    % actions
    if self.selectingEnabled
      if is1
        self.pm.updateSelectionMode(SelectionMode.Pointwise);
      end
      if is2
        self.pm.updateSelectionMode(SelectionMode.Lasso);
      end
      if isC
        self.pm.deselectAll();
      end
      switch self.pm.selectionMode
        case SelectionMode.Pointwise
          if isEsc
            self.pm.deselectLastIndex();
          end
        case SelectionMode.Lasso
          if isEsc
            self.pm.undoLassoPoint(false)
          elseif isQ
            self.pm.undoLassoPoint(true);
          elseif isEnter
            if self.pm.selectLassoedPoints();
              self.pm.undoLassoPoint(true);
            end
          end
        otherwise
          error('Bug found!');
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
  
  function refreshLasso(self)
    p = self.pm.lassoPoints; % shorthand notation
    if ~isempty(p)
      xy = self.pointerLocInCoord();
      set(self.lassoH, 'Visible', 'on', ...
        'XData', [p(:,1)' xy(1) p(1,1)],...
        'YData', [p(:,2)' xy(2) p(1,2)], 'Color', self.pm.lassoColor);
    else
      set(self.lassoH, 'Visible', 'off');
    end
  end
  
  function refresh(self)
%     if self.pm.selectionMode == SelectionMode.Lasso
%       self.refreshLasso();
% %       p = self.pm.lassoPoints; % shorthand notation
% %       if ~isempty(p)
% %         set(self.lassoH, 'Visible', 'on', 'XData', [p(:,1)' p(1,1)],...
% %           'YData', [p(:,2)' p(1,2)], 'Color', self.pm.lassoColor);
% %       end
%     else
%       set(self.lassoH, 'Visible', 'off');
%     end
    switch self.pm.selectionMode
      case SelectionMode.Pointwise
        mode = ' (pointwise select)';
        set(self.lassoH, 'Visible', 'off');
      case SelectionMode.Lasso        
        mode = ' (lasso select)';
        self.refreshLasso();
      otherwise
        error('Bug found!')
    end
    set(self.figH, 'Name', [self.pm.title mode]); 
  end
  
  function updatePlot(self)
    if self.pm.clusterCount > 8
      self.colorMap=distinguishable_colors(self.pm.clusterCount);
    else
      self.colorMap= brewermap(8,'Dark2');
    end
    if ishandle(self.figH)
      self.refresh();
      set(self.figH,'color','w');
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
        if self.pm.clusterCount==1, legend(self.figAxH,'off');
        else
          idx=unique(self.pm.cluster,'stable');
          c=unique(cmap,'rows','stable');
          s={};
          for i=1:length(idx)
              hold(self.figAxH,'on')
              h(i)=plot(self.figAxH,0,0,'o','MarkerEdgeColor',c(i,:));
              set(h(i),'Visible','off');
              if idx(i)<0
                  s{i}='scatter';
              else
                  s{i}=['Cluster ' num2str(idx(i))]; 
              end
              hold(self.figAxH,'off')
          end
          legend(h,s);
        end
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
%       colorCount = size(self.colorMap,1);
      for i = 1:size(self.pm.data,1)
        if self.pm.cluster(i)==-1
          cmap(i,:)=self.colorMap(end,:);
        else
          cmap(i,:) = self.colorMap(self.pm.cluster(i),:);
        end
      end      
    else
      cmap = self.color;    
    end
    self.cmap=cmap;
  end
  
  function updatePlotToSettings(self, settings)
    if ishandle(self.figH)
      set(self.figH, 'Position', settings);
    end
  end

  function xy = pointerLocInCoord(self)
  % Returns mouse pointer location on screen in the underlaying XY plot
  % coordinates
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

