classdef PcaTool < GuiBase & MObject
  %PCATOOL Pca tool GUI
  %   Detailed explanation goes here
  
properties (GetAccess = public, SetAccess = private)
  mainH % GUI figure
  scores % scores ScatterSelect
  loadings % loadings ScatterSelect  
  pm % Presentation model  
  settingsFile % mat-file where settings are stored
  lastSaveGenes % last save filename for genes
  lastSaveSamples % lastfilename for  samples
end


methods
  function self = PcaTool(computeObj)
  % computeObj is an object derived from PcaComputeBase
    self@MObject();
    self@GuiBase();    
    self.createScatterPlots();
    % wire presentation model signals 
    p = PcaToolPM(computeObj, self.scores.pm, self.loadings.pm);
    p.connectMe('reset', @self.refresh);
    p.connectMe('highlight_changed', @self.updateAnnotationInfo);
    p.connectMe('selection_changed', @self.updateLists);
    self.pm = p;
    self.settingsFile = fullfile(pwd, 'settings.mat');
    self.lastSaveGenes = '';
    self.lastSaveSamples = '';
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
      self.loadSettings();
      self.refresh();
      self.updateAvailableFeatures();      
      self.loadings.show();
      self.scores.show();      
      self.pm.updatePcaAxes();
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
        ind = 3 + self.pm.clusterUserInd;
      otherwise
        error('Bug found');
    end    
    str = 'Select PC axes';
    if ~isempty(self.pm.maxPcInd)
      str = [str sprintf(' (up to %d)', self.pm.maxPcInd)];
    end
    set(self.pcAxesPanelH, 'Title', str);
    str = {'k-means', 'Gaussian mixture', 'Minkowski weighted k-means'};
    str = [str self.pm.clusterUserNames];
    set(self.clusteringPopupH, 'Value', ind, 'String', str);
    self.updateAnnotationInfo();
    self.updateLists();
  end
  
  %*** Callbacks from GUI objects defined in GUIDE
  function refreshPcaButtonH_Callback(self, varargin)
    self.pm.updatePcaAxes();
  end

  function pcaxEditH_Callback(self, varargin)
    [tf, ind] = self.validatePcInput(get(self.pcaxEditH, 'String'));
    if tf
      self.pm.pcaxInd = ind;
    else
      set(self.pcaxEditH, 'String', num2str(self.pm.pcaxInd));
    end
  end
  
  function pcayEditH_Callback(self, varargin)
    [tf, ind] = self.validatePcInput(get(self.pcayEditH, 'String'));
    if tf
      self.pm.pcayInd = ind;
    else
      set(self.pcayEditH, 'String', num2str(self.pm.pcayInd));
    end    
%     ind = str2double(get(self.pcayEditH, 'String'));
%     self.pm.pcayInd = ind;
  end
  
  function clusteringPopupH_Callback(self, varargin)
    ind = get(self.clusteringPopupH, 'Value');
    if (ind <= 3)
      switch ind
        case 1
          self.pm.clusterMethod = ClusteringMethod.KMeans;
        case 2
          self.pm.clusterMethod = ClusteringMethod.Gaussian;
        case 3
          self.pm.clusterMethod = ClusteringMethod.Minkowski;
        otherwise
          error('Bug found');
      end
    else
      self.pm.clusterMethod = ClusteringMethod.User;
      self.pm.clusterUserInd = ind-3;
    end
  end
  
  function loadListsButtonH_Callback(self, varargin)
    [fname, pname, ~] = uigetfile('*.mat', 'Load lists struc');
    if fname ~= 0
      tmp = load(fullfile(pname, fname), 'cluster');
      if isempty(tmp) || ~isstruct(tmp.cluster)
        errordlg('Struct ''cluster'' not found');
      else
        self.pm.parseUserLists(tmp.cluster);
      end
    end
  end
  
  function clusterCellsButtonH_Callback(self, varargin)
    self.pm.updateCurrentClustering();    
  end
  
  function sampleListboxH_Callback(self, varargin)
    self.pm.sampleListInd = get(self.sampleListboxH, 'Value');
  end
  
  function geneListboxH_Callback(self, varargin)
    self.pm.geneListInd = get(self.geneListboxH, 'Value');
  end
  
  function clearGeneListButtonH_Callback(self, varargin)
    self.pm.deselectAllGenes();
  end
  
  function clearSampleListButtonH_Callback(self, varargin)
    self.pm.deselectAllSamples();
  end
  
  function saveSampleListButtonH_Callback(self, varargin)
    savedAs = self.openDialogAndSaveAsText('Save samples', ...
      self.lastSaveSamples, get(self.sampleListboxH, 'String'));
    if ~isempty(savedAs)
      self.lastSaveSamples = savedAs;
    end
  end
  
  function saveGeneListButtonH_Callback(self, varargin)
    savedAs = self.openDialogAndSaveAsText('Save genes', ...
      self.lastSaveGenes, get(self.geneListboxH, 'String'));
    if ~isempty(savedAs)
      self.lastSaveGenes = savedAs;
    end    
  end
  
  function deleteSampleButtonH_Callback(self, varargin)
    ind = get(self.sampleListboxH, 'Value');
    self.pm.deselectSample(ind);
  end
  
  function deleteGeneButtonH_Callback(self, varargin)
    ind = get(self.geneListboxH, 'Value');
    self.pm.deselectGene(ind)
  end
  
  function findGeneButtonH_Callback(self, varargin)
    name = strtrim(get(self.geneSymbolEditH, 'String'));
    if ~self.pm.findGene(name)
      errordlg(sprintf('Can''t find gene ''%s''', name));
    end
  end
  
  function findSampleButtonH_Callback(self, varargin)
    name = strtrim(get(self.sampleSymbolEditH, 'String'));
    if ~self.pm.findSample(name)
      errordlg(sprintf('Can''t find sample ''%s''', name));
    end
  end
  
  function addTopGenesButtonH_Callback(self, varargin)
    cutoff = str2double(get(self.cutoffEditH, 'String'))/100;
    if isnan(cutoff)
      uiwait(errordlg('Cutoff must be numeric'));
    elseif cutoff < 0 || cutoff > 1
      uiwait(errordlg('Cutoff must be between 0 and 100'));
    else
      pc = get(self.pcPopupH, 'Value');
      pos = get(self.posPopupH, 'Value');
      if pc == 1
        xOrY = 'x';
      elseif pc == 2
        xOrY = 'y';
      else
        error('Bug found');
      end
      if pos == 1
        posOrNeg = 'pos';
      elseif pos == 2
        posOrNeg = 'neg';
      else
        error('Bug found');
      end
      self.pm.selectTopGenes(cutoff, xOrY, posOrNeg);
    end
  end
  
  function refreshPcaUsingSamplesButtonH_Callback(self, varargin)
    self.pm.newPcaUsingSamples();
  end
  
  function geneSymbolEditH_Callback(self, varargin)
    % nothing to do
  end
  
  function sampleSymbolEditH_Callback(self, varargin)
    % nothing to do
  end
  
  function cutoffEditH_Callback(self, varargin)
    % nothing to do
  end
  
  function pcPopupH_Callback(self, varargin)
    % nothing to do
  end
  
  function posPopupH_Callback(self, varargin)
    % nothing to do
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
  function updateAvailableFeatures(self)
%     tags ={'saveGeneListButtonH', 'addTopGenesButtonH', ...
%       'deleteGeneButtonH', 'ontologyButtonH',...
%       'cutoffEditH', 'pc1PopupH', 'posPopupH', 'findGeneButtonH', ...
%       'geneSymbolEditH', 'deleteSampleButtonH', ...
%       'refreshPcaUsingSamplesButtonH', 'findSampleButtonH', ...
%       'sampleSymbolEditH', 'saveSampleListButtonH', 'tracePopupH', ...
%       'runTraceButtonH'};
    tags ={ 'tracePopupH', 'runTraceButtonH', 'ontologyButtonH'};
    for i = 1:length(tags)
      set(self.(tags{i}), 'Enable', 'off');
    end
  end
  
  
  function [tf, ind] = validatePcInput(self, str)
    ind = str2double(str);
    tf = false;
    if isnan(ind)    
      uiwait(errordlg('PC axis must be numeric'));
    elseif ind < 1 || ind > self.pm.maxPcInd
      uiwait(errordlg(sprintf('PC axis must be between 1 and %d', ...
        self.pm.maxPcInd)));
    else
      tf = true;
    end
  end
  
  function savedAs = openDialogAndSaveAsText(self, title, lastSave, data)
  % data is a cell array of strings
    [fname, pname, ~] = uiputfile(lastSave, title);
    if fname ~= 0
      savedAs = fullfile(pname, fname);
      fid = fopen(savedAs, 'w', 'n', 'UTF-8');
      for i=1:length(data)
        fprintf(fid, '%s\n', data{i});
      end
      fclose(fid);
    else
      savedAs = [];
    end
  end
  
  function windowAboutToClose(self, varargin)
    self.saveSettingsAndQuit();
  end
  
  function createScatterPlots(self)
  % Creates scatter plot windows, signals are fed to the presentation
  % model which handles the UI logic
    % pca scores, i.e. samples
    s = ScatterSelect();
    s.title = 'PCA scores (genes)';
    s.selectingEnabled = true;
    s.highMarker = '*';
    s.closeFcn = @self.saveSettingsAndQuit;
    self.scores = s;
    % pca loadings, i.e., genes
    s = ScatterSelect();
    s.title = 'PCA loadings (samples)';
    s.selectingEnabled = true;
    s.highMarker = '*';
    s.closeFcn = @self.saveSettingsAndQuit;
    self.loadings = s;
  end  
  
  function updateFromUiState(self)
    tags = {'deleteGeneButtonH', 'clearGeneListButtonH', ...
      'saveGeneListButtonH', 'deleteSampleButtonH', ...
      'clearSampleListButtonH', 'saveSampleListButtonH',...
      'refreshPcaUsingSamplesButtonH', 'pcaxEditH', 'pcayEditH', ....
      'refreshPcaButtonH', 'clusteringPopupH', 'clusterCellsButtonH'};
    props = {'deleteGeneEnable', 'clearGeneListEnable', ...
      'saveGeneListEnable', 'deleteSampleEnable', ...
      'clearSampleListEnable', 'saveSampleListEnable', ...
      'refreshPcaUsingSamplesEnable', 'pcxAxisEnable', 'pcyAxisEnable', ...
      'refreshPcaEnable', 'refreshPcaEnable', 'clusteringMenuEnable', ...
      'clusteringButtonEnable'};
    state = self.pm.uiState;
    for i = 1:length(tags)
      tag = tags{i};
      prop = props{i};
      if state.(prop)
        onOff = 'on';
      else
        onOff = 'off';
      end
      set(self.(tag), 'Enable', onOff);
    end
  end

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
    self.updateSampleList();
  end
  
  function updateGeneList(self)    
    ind = self.pm.geneSelIndices;
    N = length(ind);
    list = cell(N, 1);
    for i = 1:N
      list{i} = self.pm.getAnnotation('symbol_text', ind(i));
    end
    ind = self.pm.geneListInd;
    if isempty(ind)
      set(self.geneListboxH, 'String', 'Add genes to the list', ...
        'Value', 1);
    else
      set(self.geneListboxH, 'String', list, 'Value', ind);
    end
    list = cell(2,1);
    list{1} = sprintf('PC%d', self.pm.pcaxInd);
    list{2} = sprintf('PC%d', self.pm.pcayInd);
    set(self.pcPopupH, 'String', list);
    self.updateFromUiState();
  end
  
  function updateSampleList(self)
    ind = self.pm.sampleSelIndices;
    N = length(ind);
    list = cell(N, 1);
    for i = 1:N
      list{i} = self.pm.getAnnotation('id_text', ind(i));
    end
    ind = self.pm.sampleListInd;
    if isempty(ind)
      set(self.sampleListboxH, 'String', 'Add samples to the list', ...
        'Value', 1);
    else
      set(self.sampleListboxH, 'String', list, 'Value', ind);
    end
  end
  
  function saveSettings(self)
  % Save settings, e.g, figure locations, etc., when closing the tool    
    settings = struct;      
    settings.main = get(self.mainH, 'Position');
    settings.loadings = self.loadings.getSettings();
    settings.scores = self.scores.getSettings();
    settings.pm = self.pm.getSettings();
    settings.lastSaveGenes = self.lastSaveGenes;
    settings.lastSaveSamples = self.lastSaveSamples;
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
      self.lastSaveGenes = settings.lastSaveGenes;
      self.lastSaveSamples = settings.lastSaveSamples;      
    end
  end  
end
end

