classdef PcaTool < GuiBase & MObject
  %PCATOOL Pca tool GUI
  %   Detailed explanation goes here
  
properties (GetAccess = public, SetAccess = private)
  mainH % GUI figure
  scores % scores ScatterSelect
  loadings % loadings ScatterSelect  
  pm % Presentation model  
  settingsFile
end


methods
  function self = PcaTool(computeObj)
  % computeObj is an object derived from PcaComputeBase
    self@MObject();
    self@GuiBase();    
    % wire presentation model signals 
    p = PcaToolPM(computeObj);
    p.connectMe('all_changed', @self.refresh);
    p.connectMe('highlight_changed', @self.updateAnnotationInfo);
%     p.connectMe('no_highlight_changed', @()updateAnnotationInfo(self, -1));
    p.connectMe('selection_changed', @self.updateLists);
    self.pm = p;
    self.settingsFile = fullfile(pwd, 'settings.mat');
  end

  function show(self)
    if ishandle(self.mainH)
      set(self.mainH, 'Visible', 'on');
      self.registerCallbacks();
    else
      self.registerCallbacks();
      fig = pca_tool2('HandleVisiblity', 'off');
      self.mapUiControls(fig);
      set(fig, 'CloserequestFcn', @self.windowAboutToClose);
      % save callbacks
      handles = guidata(fig);
      handles.n = {self};
      guidata(fig, handles);
      self.mainH = fig;
      self.createScatterPlots();
      self.loadSettings();
      self.refresh();
      self.loadings.show();
      self.scores.show();
    end
  end

  function refresh(self)
    set(self.pcaxEditH, 'String', num2str(self.pm.pcaxInd));
    set(self.pcayEditH, 'String', num2str(self.pm.pcayInd));
    switch self.pm.clusterMethod
      case ClusteringMethod.KMeans
        ind = 1;
      case ClusteringMethod.Gaussian
        ind = 2;
      case ClusteringMethod.Minkowski
        ind = 3;
      case ClusteringMethod.User
        ind = 4;      
      otherwise
        error('Bug found');
    end    
    set(self.clusteringPopupH, 'Value', ind);
    self.refreshPcaButtonH_Callback();
    self.updateAnnotationInfo();
  end
  
  %*** Callbacks from GUI objects defined in GUIDE
  function refreshPcaButtonH_Callback(self, varargin)
    self.pm.updateCurrentPca();
    self.scores.updateData(self.pm.scoreXY, self.pm.cluster);
    self.loadings.updateData(self.pm.coefXY, self.pm.cluster);
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
        self.pm.clusterMethod = ClusteringMethod.KMeans;
      case 2
        self.pm.clusterMethod = ClusteringMethod.Gaussian;
      case 3
        self.pm.clusterMethod = ClusteringMethod.Minkowski;
      case 4
        self.pm.clusterMethod = ClusteringMethod.User;        
      otherwise
        error('Bug found');
    end
  end
  
  function clusterCellsButtonH_Callback(self, varargin)
    self.pm.updateCurrentClustering();
    self.refreshPcaButtonH_Callback();
  end
  
  %*** 
  function saveSettingsAndQuit(self)
    self.saveSettings();
    % Delete all figs
    self.loadings.closeFigure();
    self.scores.closeFigure();
    delete(self.mainH);
  end
  
end

%*** Private implementation related stuff
methods (Access = private)
  function windowAboutToClose(self, varargin)
    self.saveSettingsAndQuit();
  end
  
  function createScatterPlots(self)
  % Creates scatter plot windows, signals are fed to the presentation
  % model which handles the UI logic
    % pca scores
    s = ScatterSelect();
    s.title = 'PCA scores';
    s.selectingEnabled = true;
    s.highMarker = '*';
    s.closeFcn = @self.saveSettingsAndQuit;
    s.connectMe('is_in', @(x)registerIsIn(self.pm, 'sample', x));
    s.connectMe('highlight', @(x)highlightChanged(self.pm, 'sample', x));
    s.connectMe('selection', ...
       @(x,y)self.pm.selectionChanged('sample', x,y));
    self.scores = s;
    % pca loadings
    s = ScatterSelect();
    s.title = 'PCA loadings';
    s.selectingEnabled = true;
    s.highMarker = '*';
    s.closeFcn = @self.saveSettingsAndQuit;
    s.connectMe('is_in', @(x)registerIsIn(self.pm, 'gene', x));
    s.connectMe('highlight', @(x)highlightChanged(self.pm, 'gene', x));
    s.connectMe('selection', ...
       @(x,y)selectionChanged(self.pm, 'gene', x,y));    
    self.loadings = s;
  end

  
  %****
  
%   function clearAnnotation(self)
%     self.updateGeneAnnotations(-1);
%     self.updateCellAnnotations(-1);    
%   end
  
  function updateAnnotationInfo(self)
    self.updateGeneAnnotations(self.pm.geneHighInd);
    self.updateCellAnnotations(self.pm.sampleHighInd);
  end
  
  function updateGeneAnnotations(self, ind)
    if ~isempty(ind)
       symbolText = self.pm.getAnnotation('symbol_text', ind);
       medianText = num2str(self.pm.getAnnotation(...
         'median_number', ind));
      dispersionText = num2str(self.pm.getAnnotation(...
        'dispersion_number', ind));
      expressingText = num2str(self.pm.getAnnotation(...
        'expressing_number', ind));      
    else
      symbolText = '-';
      medianText = '-';
      dispersionText = '-';
      expressingText = '-';
    end
    set(self.symbolTextH, 'String', symbolText);
    set(self.medianTextH, 'String', medianText);
    set(self.dispersionTextH, 'String', dispersionText);
    set(self.expressingTextH, 'String', expressingText);    
  end
  
  function updateCellAnnotations(self, ind)
    titleText = 'Cell annotations';
    if ~isempty(ind)
      titleText = [titleText ' (ID ' ...
        self.pm.getAnnotation('id_text', ind) ')'];
      tagsText = num2str(self.pm.getAnnotation('tags_number', ind));
      genesText = num2str(self.pm.getAnnotation('genes_number', ind)); 
      preseqText = num2str(self.pm.getAnnotation('preseq_number', ind)); 
      simpsonText = num2str(self.pm.getAnnotation('simpson_number', ind)); 
      binomialText = num2str(self.pm.getAnnotation(...
                                              'binomial_number', ind)); 
    else
      tagsText = '-';
      genesText = '-';
      preseqText = '-';
      simpsonText = '-';
      binomialText = '-';
    end
    set(self.cellPanelH, 'Title', titleText);
    set(self.tagsTextH, 'String', tagsText);
    set(self.genesTextH, 'String', genesText);
    set(self.preseqTextH, 'String', preseqText);
    set(self.simpsonTextH, 'String', simpsonText);
    set(self.binomialTextH, 'String', binomialText);
  end
  
  function updateLists(self)
    self.updateGeneList();
  end
  
  function updateGeneList(self)    
    ind = self.pm.geneSelIndices;
    N = length(ind);
    list = cell(N, 1);
    for i = 1:N
      list{i} = self.pm.getAnnotation('symbol_text', ind(i));
    end
    set(self.geneListboxH, 'String', list);
  end
  
  function saveSettings(self)
  % Save settings, e.g, figure locations, etc., when closing the tool    
    settings = struct;      
    settings.main = get(self.mainH, 'Position');
    settings.loadings = self.loadings.getSettings();
    settings.scores = self.scores.getSettings();
    settings.pm = self.pm.getSettings();
    save(self.settingsFile, 'settings');
  end

  function loadSettings(self)
  % Load settings when the tool is opened
  % For some reason, setting the GUI window position does not work
    if exist(self.settingsFile, 'file')
      tmp = load(self.settingsFile);
      settings = tmp.settings;
      set(self.mainH, 'Position', settings.main);
      self.loadings.changeToSettings(settings.loadings);
      self.scores.changeToSettings(settings.scores);
      self.pm.changeToSettings(settings.pm);
    end
  end  
end
end

