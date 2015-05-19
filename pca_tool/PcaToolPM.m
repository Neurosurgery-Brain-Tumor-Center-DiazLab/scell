classdef PcaToolPM < handle & MObject
  %PCATOOLPM Presentation model for the PcaTool
  %   Presentation model represent the logical contents of an UI.
  %   See http://martinfowler.com/eaaDev/PresentationModel.html
  
properties
  pcaxInd % pca x index
  pcayInd % pca y index  
  clusterMethod  % current ClusteringMethod
  geneListInd
  sampleListInd  
end
  
properties (SetAccess = private, GetAccess = public)
  compute % PcaComputeBase derived object
  coefXY % current pca axes coeff values 
  scoreXY % current pca axes score values
  cluster % cluster tags
  geneHighInd % current gene highlight index, empty if no highlight
  sampleHighInd % 
  sampleIsIn = false
  geneIsIn = false  
  geneSelData  % gene selection data Mx2
  geneSelIndices % gene selection indices Mx1
  sampleSelData
  sampleSelIndices    
  uiState
end

methods
  function self = PcaToolPM(computeObj)
    self@MObject();
    computeObj.connectMe('d_changed', @self.reset);
    self.compute = computeObj;
    self.uiState = UiState();
    self.changeToSettings(self.defaultSettings);     
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
      val = -1234;      
    elseif strcmp(name, 'tags_number')
      val = sum(self.compute.d.counts(:,ind));      
    elseif strcmp(name, 'genes_number')
      val = -4444;
%       val = self.compute.d.nnz(ind,self.compute.d.counts(:,ind));
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
    self.emit('reset_changed');
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
  function registerIsIn(self, from, isIn)
    doEmit = false;
    if strcmp(from, 'sample')
      if self.sampleIsIn ~= isIn
        doEmit = true;
      end
      self.sampleIsIn = isIn;
    elseif strcmp(from, 'gene')
      if self.geneIsIn ~= isIn
        doEmit = true;
      end
      self.geneIsIn = isIn;
    else
      error('Bug found');
    end
    if ~self.sampleIsIn
      self.sampleHighInd = [];
    end
    if ~self.geneIsIn
      self.geneHighInd = [];
    end
    if doEmit
      self.emit('highlight_changed');  
    end
%     if self.sampleIsIn == false && self.geneIsIn == false
%       self.emit('no_highlight_changed');
%     end
  end  
  
  function highlightChanged(self, type, ind)
    if strcmp(type, 'gene')
      self.geneHighInd = ind;
    elseif strcmp(type, 'sample')
      self.sampleHighInd = ind;
    else
      error('Bug found');
    end    
    self.emit('highlight_changed');
  end
  
  function selectionChanged(self, type, data, indices)
    if strcmp(type, 'gene')
      self.geneSelData = data;
      self.geneSelIndices = indices;
      if isempty(indices)
        self.geneListInd = [];
      elseif isempty(self.geneListInd) && ~isempty(indices)
        self.geneListInd = 1;
      elseif self.geneListInd > length(indices)
        self.geneListInd = length(indices);
      end
      if isempty(self.geneListInd)
        self.uiState.updateGenesSelected(false);
      else
        self.uiState.updateGenesSelected(true);
      end
    elseif strcmp(type, 'sample')
      self.sampleSelData = data;
      self.sampleSelIndices = indices;
      if isempty(indices)
        self.sampleListInd = [];
      elseif isempty(self.sampleListInd) && ~isempty(indices)
        self.sampleListInd = 1;        
      elseif self.sampleListInd > length(indices)
        self.sampleListInd = length(indices);
      end
      if isempty(self.sampleListInd)
        self.uiState.updateSamplesSelected(false);
      else
        self.uiState.updateSamplesSelected(true);
      end
    else
      error('Bug found');
    end
    self.emit('selection_changed')
  end
  
end

methods (Access = private)

end
  
end

