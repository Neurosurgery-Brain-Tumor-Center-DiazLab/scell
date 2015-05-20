classdef PcaToolPM < MObject
  %PCATOOLPM Presentation model for the PcaTool
  %   Presentation model represent the logical contents of an UI.
  %   See http://martinfowler.com/eaaDev/PresentationModel.html
  
properties
  pcaxInd % pca x index
  pcayInd % pca y index  
  clusterMethod  % current ClusteringMethod
  geneListInd % current selected index in gene list
  sampleListInd % current selected index in sample list 
end
  
properties (SetAccess = private, GetAccess = public)
  compute % PcaComputeBase derived object
  coefXY % current pca axes coeff values 
  scoreXY % current pca axes score values
  cluster % cluster tags
  samplePm % presentation model for samples, i.e., scores window 
  genePm % presentation model for genes, i.e., loadings window
 
  uiState
end

properties (Dependent)
  % dep
  geneHighInd % current gene highlight index, empty if no highlight
  sampleHighInd % 
%   sampleIsIn = false
%   geneIsIn = false  
  geneSelData  % gene selection data Mx2
  geneSelIndices % gene selection indices Mx1
  sampleSelData
  sampleSelIndices     
end

methods
  function self = PcaToolPM(computeObj, samplePm, genePm)
    self@MObject();
    computeObj.connectMe('d_changed', @self.reset);
    % just forward (re-emit) selection signals
    samplePm.connectMe('highlight_changed',... 
        @(x)emit(self, 'highlight_changed'));
    genePm.connectMe('highlight_changed',... 
        @(x)emit(self, 'highlight_changed'));
    samplePm.connectMe('selection_changed', @self.sampleSelectionChanged);
    genePm.connectMe('selection_changed', @self.geneSelectionChanged);     
    self.compute = computeObj;
    self.samplePm = samplePm;
    self.genePm = genePm;
    self.uiState = UiState();
    self.changeToSettings(self.defaultSettings);     
  end

  function val = get.geneHighInd(self)
    val = self.genePm.currentIndex;
  end
  
  function val = get.sampleHighInd(self)
    val = self.samplePm.currentIndex;
  end
  
  function val = get.geneSelData(self)
    val = self.genePm.selData;
  end
  
  function val = get.geneSelIndices(self)
    val = self.genePm.selIndices;
  end
  
  function val = get.sampleSelData(self)
    val = self.samplePm.selData;
  end
  
  function val = get.sampleSelIndices(self)
    val = self.samplePm.selIndices;
  end
  
  function val = defaultSettings(self)
    val.pcaxInd = 1;
    val.pcayInd = 2;
    val.clusterMethod = ClusteringMethod.KMeans;
  end
  
  function val = getSettings(self)
    val.pcaxInd = self.pcaxInd;
    val.pcayInd = self.pcayInd;
    val.clusterMethod = self.clusterMethod;
  end
  
  function val = getAnnotation(self, name, ind)
    if strcmp(name, 'symbol_text')
      val = self.compute.d.gsymb{ind};
    elseif strcmp(name, 'id_text')
      val = self.compute.d.slbls{ind};
    elseif strcmp(name, 'median_number')
      val = median(self.compute.d.cpm(ind,:));
    elseif strcmp(name, 'dispersion_number')
      val = self.compute.d.iod(ind);
    elseif strcmp(name, 'expressing_number')
      val =  nnz(self.compute.d.counts(ind,:));
    elseif strcmp(name, 'tags_number')
      val = sum(self.compute.d.counts(:,ind));      
    elseif strcmp(name, 'genes_number')
      val = nnz(self.compute.d.counts(:,ind));
    elseif strcmp(name, 'preseq_number')
      val = self.compute.d.preseq(ind);
    elseif strcmp(name, 'simpson_number')
      val = self.compute.d.simpson(ind);
    elseif strcmp(name, 'binomial_number')
      val = self.compute.d.lorenz(ind);
    else
      error('Unknown name %s', name);
    end
  end
  
  function reset(self)
    self.uiState = UiState();
    self.emit('reset');
  end 
  
  function changeToSettings(self, settings)
    self.pcaxInd = settings.pcaxInd;
    self.pcayInd = settings.pcayInd;
    self.clusterMethod = settings.clusterMethod;
  end
  
  function updateCurrentPca(self)
    self.compute.computePca();
    if ~isempty(self.compute.coeff)
      self.coefXY = [self.compute.coeff(:,self.pcaxInd) ...
                      self.compute.coeff(:,self.pcayInd)];
      self.scoreXY = [self.compute.score(:,self.pcaxInd) ...
                      self.compute.score(:,self.pcayInd)];                  
    end
  end
  
  function updateCurrentClustering(self)
    if ~isempty(self.coefXY)
      self.cluster = randi(double(self.clusterMethod), size(self.coefXY,1));
    end
  end 
  
  %**** Handlers for changes in the UI 
  function sampleSelectionChanged(self)
    % move the current selected index in the list, if needed
    indices = self.sampleSelIndices;
    if isempty(indices)
      self.sampleListInd = [];
    elseif isempty(self.sampleListInd) && ~isempty(indices)
      self.sampleListInd = 1;        
    elseif self.sampleListInd > length(indices)
      self.sampleListInd = length(indices);
    end
    % update UI state
    if isempty(self.sampleListInd)
      self.uiState.updateSamplesSelected(false);
    else
      self.uiState.updateSamplesSelected(true);
    end    
    self.emit('selection_changed');
  end
  
  function geneSelectionChanged(self)
    % move the current selected index in the list, if needed
    indices = self.geneSelIndices;
    if isempty(indices)
      self.geneListInd = [];
    elseif isempty(self.geneListInd) && ~isempty(indices)
      self.geneListInd = 1;        
    elseif self.geneListInd > length(indices)
      self.geneListInd = length(indices);
    end
    % update UI states
    if isempty(self.geneListInd)
      self.uiState.updateGenesSelected(false);
    else
      self.uiState.updateGenesSelected(true);
    end    
    self.emit('selection_changed');
  end
  
  function deselectSample(self, index)
    self.samplePm.deselectIndex(index);
  end
  
  function deselectGene(self, index)
    self.genePm.deselectIndex(index);
  end
  
  function deselectAllSamples(self)
    self.samplePm.deselectAll();
  end
  
  function deselectAllGenes(self)
    self.genePm.deselectAll();
  end
  
%   function registerIsIn(self, from, isIn)
%     doEmit = false;
%     if strcmp(from, 'sample')
%       if self.sampleIsIn ~= isIn
%         doEmit = true;
%       end
%       self.sampleIsIn = isIn;
%     elseif strcmp(from, 'gene')
%       if self.geneIsIn ~= isIn
%         doEmit = true;
%       end
%       self.geneIsIn = isIn;
%     else
%       error('Bug found');
%     end
%     if ~self.sampleIsIn
%       self.sampleHighInd = [];
%     end
%     if ~self.geneIsIn
%       self.geneHighInd = [];
%     end
%     if doEmit
%       self.emit('highlight_changed');  
%     end
%     if self.sampleIsIn == false && self.geneIsIn == false
%       self.emit('no_highlight_changed');
%     end
%   end  
  
%   function highlightChanged(self, type, ind)
%     if strcmp(type, 'gene')
%       self.geneHighInd = ind;
%     elseif strcmp(type, 'sample')
%       self.sampleHighInd = ind;
%     else
%       error('Bug found');
%     end    
%     self.emit('highlight_changed');
%   end
%   

%   function selectionChanged(self, type, data, indices)
%     if strcmp(type, 'gene')
%       self.geneSelData = data;
%       self.geneSelIndices = indices;
%       if isempty(indices)
%         self.geneListInd = [];
%       elseif isempty(self.geneListInd) && ~isempty(indices)
%         self.geneListInd = 1;
%       elseif self.geneListInd > length(indices)
%         self.geneListInd = length(indices);
%       end
%       if isempty(self.geneListInd)
%         self.uiState.updateGenesSelected(false);
%       else
%         self.uiState.updateGenesSelected(true);
%       end
%     elseif strcmp(type, 'sample')
%       self.sampleSelData = data;
%       self.sampleSelIndices = indices;
%       if isempty(indices)
%         self.sampleListInd = [];
%       elseif isempty(self.sampleListInd) && ~isempty(indices)
%         self.sampleListInd = 1;        
%       elseif self.sampleListInd > length(indices)
%         self.sampleListInd = length(indices);
%       end
%       if isempty(self.sampleListInd)
%         self.uiState.updateSamplesSelected(false);
%       else
%         self.uiState.updateSamplesSelected(true);
%       end
%     else
%       error('Bug found');
%     end
%     self.emit('selection_changed')
%   end
  
end

methods (Access = private)

end
  
end

