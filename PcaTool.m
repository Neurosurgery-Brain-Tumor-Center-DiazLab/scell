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
    self.settingsFile = fullfile(ctfroot, 'settings.mat');
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
    set(self.plotPopupH, 'String', {'contour','surface'});
    set(self.plotMethodPopupH, 'String', {'Lowess','Loess','Linear','Cubic','Bi-harmonic','Thin-plate'});
    switch self.pm.clusterMethod
      case ClusteringMethod.KMeans
        ind = 1;
      case ClusteringMethod.Gaussian
        ind = 2;
      case ClusteringMethod.Minkowski
        ind = 3;
      case ClusteringMethod.DBSCAN
        ind = 4;
      case ClusteringMethod.User
        ind = 5;
      otherwise
        error('Bug found');
    end    
    str = 'Select PC axes';
    if ~isempty(self.pm.maxPcInd)
      str = [str sprintf(' (up to %d)', self.pm.maxPcInd)];
    end
    set(self.pcAxesPanelH, 'Title', str);
    str = {'k-means', 'Gaussian mixture', 'Minkowski weighted k-means', 'DBSCAN','User'};
    str = [str self.pm.clusterUserNames];
    set(self.clusteringPopupH, 'Value', ind, 'String', str);
    idx=unique(self.pm.cluster);
    if isempty(idx)||length(idx)==1
        set(self.setRootPopupH,'Value',1,'String','no clusters');
    else
        s={};k=1;
        for i=1:length(idx)
            if idx(i)>0
                s{k}=['Cluster ' num2str(idx(i))];
                k=k+1;
            end
        end
        set(self.setRootPopupH,'Value',1,'String',s);
    end
    self.updateAnnotationInfo();
    self.updateLists();
  end
  
  %*** Callbacks from GUI objects defined in GUIDE
  function refreshPcaButtonH_Callback(self, varargin)
    self.pm.updateCurrentClustering([],[]);
    self.pm.updatePcaAxes();
  end
  
  function varimaxRotateButtonH_Callback(self, varargin)
    cnew=self.pm.compute.coeff;
    snew=self.pm.compute.score;
    i1=self.pm.pcaxInd; i2=self.pm.pcayInd;
    [c,R]=rotatefactors(cnew(:,[i1,i2]),'Method','varimax','Normalize','off');
    cnew(:,[i1,i2])=c;
    snew(:,[i1,i2])=snew(:,[i1,i2])*R;
    self.pm.compute.changeCoeff(cnew);
    self.pm.compute.changeScore(snew);
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
  end
  
  function clusteringPopupH_Callback(self, varargin)
    ind = get(self.clusteringPopupH, 'Value');
    if (ind <= 4)
      switch ind
        case 1
          self.pm.clusterMethod = ClusteringMethod.KMeans;
        case 2
          self.pm.clusterMethod = ClusteringMethod.Gaussian;
        case 3
          self.pm.clusterMethod = ClusteringMethod.Minkowski;
        case 4
          self.pm.clusterMethod = ClusteringMethod.DBSCAN
        case 5
          self.pm.clusterMethod = ClusteringMethod.User;
        otherwise
          error('Bug found');
      end
    else
      self.pm.clusterMethod = ClusteringMethod.User;
      self.pm.clusterUserInd = ind-5;
    end
  end
  
%   function loadListsButtonH_Callback(self, varargin)
%     [fname, pname, ~] = uigetfile('*.mat', 'Load lists struc');
%     if fname ~= 0
%       tmp = load(fullfile(pname, fname), 'cluster');
%       if isempty(tmp) || ~isstruct(tmp.cluster)
%         errordlg('Struct ''cluster'' not found');
%       else
%         self.pm.parseUserLists(tmp.cluster);
%       end
%     end
%   end
  
  function clusterCellsButtonH_Callback(self, varargin)
    if ~isempty(self.pm.coefXY)
      switch self.pm.clusterMethod
          case ClusteringMethod.KMeans
              km_param=choose_kmeans_params('String','Set k-means parameters...');
              if isempty(gcp('nocreate')), parpool; end
              opt=statset('UseParallel','always');
              h=waitbar(0.5,'Clustering...');
              [U,C]=kmeans(self.pm.scoreXY,km_param.num_clust,'Replicates',km_param.num_reps,'Distance',km_param.dist,'EmptyAction','drop','Options',opt);
          case ClusteringMethod.Gaussian
              gm_param=choose_gauss_params('String','Set Gaussian MM parameters...');
              if isempty(gcp('nocreate')), parpool; end
              h=waitbar(0.5,'Clustering...');
              opt=statset('UseParallel',true,'MaxIter',1000);
              gm=fitgmdist(self.pm.scoreXY,gm_param.num_clust,'Replicates',gm_param.num_reps,'Options',opt);
              U=cluster(gm,self.pm.scoreXY);
              for i=1:max(U), C(i,:)=mean(self.pm.scoreXY(find(U==i),:)); end
          case ClusteringMethod.Minkowski
              lnk=eye(size(self.pm.scoreXY,1));
              mk_param=choose_mink_params('String','Set Minkowski-weighted k-means parameters...');
              if isempty(gcp('nocreate')), parpool; end
              h=waitbar(0.5,'Clustering...');
              [U,~,C,~,~]=SubWkMeans(self.pm.scoreXY,mk_param.num_clust,mk_param.weight_exp,false,false,mk_param.mink_exp,lnk);
          case ClusteringMethod.DBSCAN
              dbscan_param=choose_dbscan_params('String','Set DBSCAN parameters...');
              if isempty(gcp('nocreate')), parpool; end
              h=waitbar(0.5,'Clustering...');
              [U,t]=dbscan(self.pm.scoreXY,dbscan_param.min_size,[]);
              for i=1:max(U), C(i,:)=mean(self.pm.scoreXY(find(U==i),:)); end
          case ClusteringMethod.User
              d=self.pm.compute.d;
              [fname pname]=uigetfile('*.*','Select clustering definition file...');
              try
                  f=fopen(fullfile(pname,fname));
                  D=textscan(f,'%s%n','Headerlines',1);
              catch me
                  alert('String','Error reading file');
                  return;
              end
              U=(-1)*ones(size(self.pm.scoreXY,1),1);
              for i=1:length(D{1})
                  t=find(strcmp(D{1}{i},d.slbls));
                  if ~isempty(t), U(t)=D{2}(i); end
              end
              for i=1:max(U), C(i,:)=mean(self.pm.scoreXY(find(U==i),:)); end
      end
    end  
    delete(h);
    self.pm.updateCurrentClustering(U,C);  
    self.pm.clustCmap=self.scores.cmap;
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
  
  function runTraceButtonH_Callback(self, varargin)
    d=self.pm.compute.d;
    %find samples currently displayed
    idx=[];
    for i=1:self.pm.compute.sampleCount
      t=min(find(strcmp(self.pm.getAnnotation('id_text',i),d.slbls_full)));
      idx=[idx;t];
    end
    %regress/fit a surface to gene expression for user's gene and plot
    s=get(self.genePlotEditH,'String');
    t=find(strcmp(s,d.gsymb_full));
    if isempty(t), alert('String','Gene not found!!!'); return; end
    z=max(0,log2(d.cpm_full(t,idx)'));
    fstr='lowess';
    switch get(self.plotMethodPopupH,'Value')
        case 1
            fstr='lowess';
        case 2
            fstr='loess';
        case 3
            fstr='linearinterp';
        case 4
            fstr='cubicinterp';
        case 5
            fstr='biharmonicinterp';
        case 6
            fstr='thinplateinterp';
    end
    surffit=fit(self.scores.pm.data,z,fstr,'normalize','off');
    figure;
    if get(self.plotPopupH,'Value')==1
      plot(surffit,self.scores.pm.data,z,'Style','contour')
      title(s)
      h=colorbar;
      title(h,'$\log_2$CPM','Interpreter','latex')
    else
      plot(surffit,self.scores.pm.data,z);
      title(s)
      h=colorbar
      title(h,'$\log_2$CPM','Interpreter','latex')
    end
    if ~isempty(self.pm.Tr)
      Tr=self.pm.Tr;pred=self.pm.pred
      C=self.pm.clusterCtrs;
      for i=1:length(pred)
        if pred(i)==0, continue; end
        path=graphpred2path(pred,i);
        lbls={};tk=1:length(path);
        for i=1:length(path), lbls{i}=['C' num2str(path(i))]; end
        f=figure;ax=gca;
        set(f,'color','w');
        xx=[];yy=[];t=[];
        for j=1:length(path)-1
          xq=linspace(C(path(j),1),C(path(j+1),1),10);
          xq=xq(1:end-1);
          t=[t,linspace(j,j+1,10)];
          t=t(1:end-1);
          xx=[xx,xq];
          vq = interp1([C(path(j),1),C(path(j+1),1)],[C(path(j),2),C(path(j+1),2)],xq);
          yy=[yy,vq];
        end
        zz=surffit(xx,yy);
        plot(ax,t,zz);
        title(s)
        xlabel('Cluster'); ylabel('$\log_2$CPM','Interpreter','latex');
        set(ax,'XTick',tk,'XTickLabel',lbls,'FontSize',16);
      end
    end
        
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
  
  function exportClustButtonH_Callback(self, varargin)
    [fname pname]=uiputfile('*.txt','Save data as...','clusters.txt');
    try
      f=fopen(fullfile(pname,fname),'w');
    catch me
      alert('String',['Error opening ' fullfile(pname,fname)]);
      return;
    end
    fprintf(f,'Sample_ID\tCluster_ID\n');
    slbls=self.pm.compute.d.slbls;
    for i=1:length(slbls)
      fprintf(f,'%s\t',slbls{i});
      fprintf(f,'%i\n',self.pm.cluster(i));
    end
    fclose(f);  
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
  
  function fitTreeButtonH_Callback(self, varargin)
    if isempty(self.pm.cluster)||length(self.pm.cluster)==1
      alert('String','Cluster the data first');
      return;
    end
    cmap = self.pm.clustCmap;
    s={};
    for i=1:length(self.pm.cluster)
      if self.pm.cluster(i)<0
        s{i}='scatter';
      else
        s{i}=['Cluster ' num2str(self.pm.cluster(i))]; 
      end
    end
    f=figure;
    set(f,'color','w');
    gscatter(self.pm.scoreXY(:,1),self.pm.scoreXY(:,2),s',unique(cmap,'rows','stable'),'.',20);
    ax=gca;
    hold(ax,'on')
    cD=squareform(pdist(self.pm.clusterCtrs));
    [Tr,pred]=graphminspantree(sparse(cD),get(self.setRootPopupH,'Value'));
    self.pm.Tr=Tr; self.pm.pred=pred;
    C=self.pm.clusterCtrs;
    for i=1:size(C,1)
      for j=1:size(C,1)
        if Tr(i,j)~=0
            plot([C(i,1),C(j,1)],[C(i,2),C(j,2)],'LineWidth',2)
        end
      end
    end
  end
  
  function ontologyButtonH_Callback(self, varargin)
    data=get(self.geneListboxH, 'String');
    if isempty(data)
        alert('String','Select genes first')
        return;
    end
    smp.species='hsa';
    smp.gsymb=data;
    david_tool(smp);
  end
  
  function addTopGenesButtonH_Callback(self, varargin)
    if isempty(get(self.cutoffEditH, 'String')), return;end
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
    self.pm.updateCurrentClustering([],[]);
    self.pm.newPcaUsingSamples();
  end
  
   function genePlotEditH_Callback(self, varargin)
    % nothing to do
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
  
  function plotPopupH_Callback(self, varargin)
    % nothing to do
  end
  
  function posPopupH_Callback(self, varargin)
    % nothing to do
  end
  
  function setRootPopupH_Callback(self, varargin)
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
    tags ={};
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
    s.title = 'Sample scores';
    s.selectingEnabled = true;
    s.highMarker = '*';
    s.closeFcn = @self.saveSettingsAndQuit;
    self.scores = s;
    % pca loadings, i.e., genes
    s = ScatterSelect();
    s.title = 'Gene loadings';
    s.selectingEnabled = true;
    s.highMarker = '*';
    s.closeFcn = @self.saveSettingsAndQuit;
    self.loadings = s;
  end  
  
  function updateFromUiState(self)
    % Select PC axes
    tags = {'pcaxEditH', 'pcayEditH'};
    props = {'pcxAxisEnable', 'pcyAxisEnable'};
    % Clustering
    tags = [tags {'clusteringPopupH', ...
      'clusterCellsButtonH'}];
    props = [props {'clusteringMenuEnable', 'customListEnable', ...
      'clusteringButtonEnable'}];
    % Working gene list
    tags = [tags {'deleteGeneButtonH', 'clearGeneListButtonH', ...
      'findGeneButtonH', 'geneSymbolEditH', 'saveGeneListButtonH'}];
    props = [props {'deleteGeneEnable', 'clearGeneListEnable', ...
      'findGeneEnable', 'geneSymbolInputEnable', 'saveGeneListEnable'}];
    % Working sample list
    tags = [tags {'deleteSampleButtonH', 'clearSampleListButtonH', ...
      'findSampleButtonH', 'sampleSymbolEditH', ...
      'saveSampleListButtonH', 'refreshPcaUsingSamplesButtonH'}];
    props = [props {'deleteSampleEnable', 'clearSampleListEnable', ...
      'findSampleEnable', 'sampleInputEnable', 'saveSampleListEnable', ...
      'refreshPcaUsingSamplesEnable'}];
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
      ctypeText = self.pm.getAnnotation('ctype_number', ind);
      genesText = num2str(self.pm.getAnnotation('genes_number', ind)); 
      preseqText = num2str(self.pm.getAnnotation('preseq_number', ind)); 
      simpsonText = num2str(self.pm.getAnnotation('simpson_number', ind)); 
      binomialText = num2str(self.pm.getAnnotation(...
                                              'binomial_number', ind)); 
    else
      ctypeText = '-';
      genesText = '-';
      preseqText = '-';
      simpsonText = '-';
      binomialText = '-';
    end
    set(self.cellPanelH, 'Title', titleText);
    set(self.ctypeTextH, 'String', ctypeText);
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

