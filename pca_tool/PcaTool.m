classdef PcaTool < GuiBase & MObject
  %PCATOOL Pca tool GUI
  %   Detailed explanation goes here
  
properties (GetAccess = public, SetAccess = private)
  mainH % GUI figure
  scores % scores ScatterSelect
  loadings % loadings ScatterSelect  
  pm % Presentation model
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
    else
      self.registerCallbacks();
      fig = pca_tool2();
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
    disp('ok')
  end
  
  %*** Callbacks from GUI objects defined in GUIDE
  function refreshPcaButtonH_Callback(self, varargin)
    pm.updateCurrentPca();
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
    self.scores = s;
    % pca loadings
    s = ScatterSelect();
    s.title = 'PCA loadings';
    s.closeFcn = @self.saveSettingsAndQuit;
    self.loadings = s;
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

