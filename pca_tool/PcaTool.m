classdef PcaTool < GuiBase & MObject
  %PCATOOL Pca tool GUI
  %   Detailed explanation goes here
  
properties (GetAccess = public, SetAccess = private)
  mainH % GUI figure
  scores % scores ScatterSelect
  loadings % loadings ScatterSelect  
  pm % Presentation model
  scoresIsIn = false
  loadingsIsIn = false
end

properties (Constant)
  settingsFile = 'settings.mat'
end

methods
  function self = PcaTool(computeObj)
  % computeObj is an object derived from PcaComputeBase
    self@MObject();
    self@GuiBase();    
    p = PcaToolPM(computeObj);
    p.connectMe('all_changed', @self.refresh);
    self.pm = p;
  end

  function show(self)
    if ishandle(self.mainH)
      set(self.mainH, 'Visible', 'on');
      self.registerCallbacks();
    else
      self.registerCallbacks();
      fig = pca_tool2();
      self.mapUiControls(fig);
      set(fig, 'CloserequestFcn', @self.saveSettingsAndQuit);
      % save callbacks
      handles = guidata(fig);
      handles.n = {self};
      guidata(fig, handles);
      self.mainH = fig;
      self.createScatterPlots();
      self.loadSettings();
      self.loadings.show();
      self.scores.show();
    end
  end

  function refresh(self)
    self.refreshPcaButtonH_Callback();    
  end
  
  %*** Callbacks from GUI objects defined in GUIDE
  function refreshPcaButtonH_Callback(self, varargin)
    self.pm.updateCurrentPca();
    self.scores.data = self.pm.coefXY;
    self.loadings.data = self.pm.scoreXY;
  end

  function pcaxEditH_Callback(self, varargin)
    ind = str2double(get(self.pcaxEditH, 'String'));
    self.pm.pcaxInd = ind;
  end
  
  function pcayEditH_Callback(self, varargin)
    ind = str2double(get(self.pcayEditH, 'String'));
    self.pm.pcayInd = ind;
  end
  
  function clusteringPopupH_Callback(self, varargin)
    ind = get(self.clusteringPopupH, 'Value');
    switch ind
      case 1
        self.pm.clustering = ClusteringMethod.KMeans;
      case 2
        self.pm.clustering = ClusteringMethod.Gaussian;
      case 3
        self.pm.clustering = ClusteringMethod.Minkowski;
      case 4
        self.pm.clustering = ClusteringMethod.User;        
      otherwise
        error('Bug found');
    end
  end
  
  function clusterCellsButtonH_Callback(self, varargin)
    
  end
  
  %*** Other callbacks
  function saveSettingsAndQuit(self, varargin)
    self.saveSettings();
    % Delete all figs
    self.loadings.closeFigure();
    self.scores.closeFigure();
    delete(self.mainH);
  end
end

%*** Private implementation related stuff
methods (Access = private)
  function createScatterPlots(self)
  % Creates scatter plot windows
    % pca scores
    s = ScatterSelect();
    s.title = 'PCA scores';
    s.closeFcn = @self.saveSettingsAndQuit;
    s.connectMe('is_in', @(x)self.registerIsIn('scores', x));
    s.connectMe('highlight', @self.updateAnnotationInfo);
    self.scores = s;
    % pca loadings
    s = ScatterSelect();
    s.title = 'PCA loadings';
    s.closeFcn = @self.saveSettingsAndQuit;
    s.connectMe('is_in', @(x)self.registerIsIn('loadings', x));
    s.connectMe('highlight', @self.updateAnnotationInfo);
    self.loadings = s;
  end

  function registerIsIn(self, from, isIn)
    if strcmp(from, 'scores')
      self.scoresIsIn = isIn;
    else
      self.loadingsIsIn = isIn;
    end
    if self.scoresIsIn == false && self.loadingsIsIn == false
      self.updateAnnotationInfo(-1);
    end
  end
  
  function updateAnnotationInfo(self, ind)
    self.updateGeneAnnotations(ind);
    self.updateCellAnnotations(ind);
  end
  
  function updateGeneAnnotations(self, ind)
    if ind > 0
       symbolText = self.pm.getAnnotation('symbol_text', ind);
%        medianText = self.pm.getAnnotation('median_value', ind);
      medianText = 'todo';
      dispersionText = num2str(self.pm.getAnnotation(...
        'dispersion_number', ind));
    else
      symbolText = '-';
      medianText = '-';
      dispersionText = '-';
    end
    set(self.symbolTextH, 'String', symbolText);
    set(self.medianTextH, 'String', medianText);
    set(self.dispersionTextH, 'String', dispersionText);
  end
  
  function updateCellAnnotations(self, ind)
    titleText = 'Cell annotations';
    if ind > 0
      titleText = [titleText ' (ID ' num2str(...
        self.pm.getAnnotation('id_number', ind)) ')'];
    else
      
    end
    set(self.cellPanelH, 'Title', titleText);
  end
  
  function saveSettings(self)
  % Save settings, e.g, figure locations, etc., when closing the tool    
    settings = struct;      
    settings.main = get(self.mainH, 'Position');
    settings.loadings = self.loadings.getSettings();
    settings.scores = self.scores.getSettings();
    settings.pm = self.pm.getSettings();
    save(PcaTool.settingsFile, 'settings');
  end

  function loadSettings(self)
  % Load settings when the tool is opened
  % For some reason, setting the GUI window position does not work
    if exist(PcaTool.settingsFile, 'file')
      tmp = load(PcaTool.settingsFile);
      settings = tmp.settings;
      set(self.mainH, 'Position', settings.main);
      self.loadings.changeToSettings(settings.loadings);
      self.scores.changeToSettings(settings.scores);
      self.pm.changeToSettings(settings.pm);
    end
  end  
end
end

