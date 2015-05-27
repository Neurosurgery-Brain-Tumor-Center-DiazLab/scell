classdef PcaToolPM < MObject
  %PCATOOLPM Presentation model for the PcaTool
  %   Presentation model represent the logical contents of an UI.
  %   See http://martinfowler.com/eaaDev/PresentationModel.html
  
properties
  pcaxInd % pca x index, i.e., left side PC in the GUI
  pcayInd % pca y index, i.e., right side PC in the GUI 
  clusterMethod  % current ClusteringMethod
  geneListInd % current selected index in gene list
  sampleListInd % current selected index in sample list 
end

properties (SetAccess = private, GetAccess = public)
  compute % PcaComputeBase derived object
  coefXY % current pca axes coeff values, i.e., genes
  scoreXY % current pca axes score values, i.e., samples
  cluster % cluster tags
  samplePm % presentation model for samples, i.e., scores window 
  genePm % presentation model for genes, i.e., loadings window
  nBins = 20 % number of bins in gene histogram
  geneCumProbPosX % Nx2 cumulative probability function for positive genes
                 % each row is [loading cumulative_prob]
  geneCumProbNegX % as above, but for negative value, each row is
                 % [abs(loading) cumulative_prob]
  geneCumProbPosY
  geneCumProbNegY
  uiState
  enableSelectionChanged = true % helper boolean
end

properties (Dependent)
  maxPcInd % maximum PC index, i.e, number of colums in the coeff
  geneHighInd % current gene highlight index, empty if no highlight
  sampleHighInd % 
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

  function val = get.maxPcInd(self)
    val = [];
    if ~isempty(self.compute.coeff)
      val = size(self.compute.coeff,2);
    end
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
    self.pcaxInd = 1;
    self.pcayInd = 2;    
    self.updateCurrentPca();
  end 
  
  function changeToSettings(self, settings)
    self.pcaxInd = settings.pcaxInd;
    self.pcayInd = settings.pcayInd;
    self.clusterMethod = settings.clusterMethod;
  end
  
  function updateCurrentPca(self)
    self.compute.computePca();
    if ~isempty(self.compute.coeff)
      self.uiState.updateHasData(true);
      self.coefXY = [self.compute.coeff(:,self.pcaxInd) ...
                      self.compute.coeff(:,self.pcayInd)];
      self.scoreXY = [self.compute.score(:,self.pcaxInd) ...
                      self.compute.score(:,self.pcayInd)];   
      self.samplePm.updateData(self.scoreXY, self.cluster);
      self.genePm.updateData(self.coefXY, self.cluster);
            
      % calculate cumulative probability function from positive genes X
      ind = self.coefXY(:,self.pcaxInd) >= 0;
      [n, edges] = histcounts(self.coefXY(ind,self.pcaxInd), ...
        self.nBins);
      ncum = cumsum(n);
      binCenters = edges(1:end-1) + diff(edges);
      self.geneCumProbPosX = [binCenters(:) ncum(:)/ncum(end)];
      % the same for negatives X
      ind = self.coefXY(:,self.pcaxInd) < 0;
      [n, edges] = histcounts(abs(self.coefXY(ind,self.pcaxInd)), ...
        self.nBins);
      ncum = cumsum(n);
      binCenters = edges(1:end-1) + diff(edges);
      self.geneCumProbNegX = [binCenters(:) ncum(:)/ncum(end)];      
      % posive Y
      ind = self.coefXY(:,self.pcayInd) >= 0;
      [n, edges] = histcounts(self.coefXY(ind,self.pcayInd), ...
        self.nBins);
      ncum = cumsum(n);
      binCenters = edges(1:end-1) + diff(edges);
      self.geneCumProbPosY = [binCenters(:) ncum(:)/ncum(end)];
      % negative Y
      ind = self.coefXY(:,self.pcayInd) < 0;
      [n, edges] = histcounts(abs(self.coefXY(ind,self.pcayInd)), ...
        self.nBins);
      ncum = cumsum(n);
      binCenters = edges(1:end-1) + diff(edges);
      self.geneCumProbNegY = [binCenters(:) ncum(:)/ncum(end)];         
      self.emit('reset');
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
    if self.enableSelectionChanged
      self.emit('selection_changed')
    end;
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
    if self.enableSelectionChanged
      self.emit('selection_changed');
    end
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
  
  function tf = findGene(self, name)
    tf = false;
    % loop all gene symbols, break on first match
    for i = 1:self.compute.geneCount
      if strcmp(name, self.getAnnotation('symbol_text', i))
        tf = true;
        break;
      end
    end
    if tf
      % If the index was not already in the selection list, add it
      % there. This which causes selection_changed signal, but we don't
      % want to re-emit it yet, so it's temporarily disabled.
      % After the index addtion, find a new list index and emit 
      % selection change.
      listInd = find(i == self.geneSelIndices);
      if isempty(listInd)
        self.enableSelectionChanged = false;
        self.genePm.selectIndex(i);
        self.enableSelectionChanged = true;
        listInd = find(i == self.geneSelIndices);
      end
      self.geneListInd = listInd;
      self.emit('selection_changed');
    end    
  end
  
  function tf = findSample(self, name)
    % see findGene for comments
    tf = false;
    for i = 1:self.compute.sampleCount
      if strcmp(name, self.getAnnotation('id_text', i))
        tf = true;
        break;
      end
    end
    if tf
      listInd = find(i == self.sampleSelIndices);
      if isempty(listInd)
        self.enableSelectionChanged = false;
        self.samplePm.selectIndex(i);
        self.enableSelectionChanged = true;
        listInd = find(i == self.sampleSelIndices);
      end
      self.sampleListInd = listInd;
      self.emit('selection_changed');
    end   
  end
  
  function limit = solveLoadingLimit(self, geneCumProb, prob)
    if prob <= eps 
      limit = geneCumProb(1,1);
    elseif prob >= 1-eps
      limit = geneCumProb(end,1);
    else
      % function whose root is to be found
      f = @(x)(interp1(geneCumProb(:,1), geneCumProb(:,2),x) - prob);
      % solution is to be found within this limit
      init = [geneCumProb(1,1) geneCumProb(end,1)];
      % find the root
      limit = fzero(f, init);
    end
  end
  
  function selectTopGenes(self, cutoff, xOrY, posOrNeg)
  % 0<= cutoff <= 1
    if ~ismember(xOrY, {'x', 'y'}) || ~ismember(posOrNeg, {'pos', 'neg'})
      error('Bug found');
    end

    % select cum prob distribution
    if strcmp(posOrNeg, 'pos') && strcmp(xOrY, 'x')
      geneCumProb = self.geneCumProbPosX;     
    elseif strcmp(posOrNeg, 'pos') && strcmp(xOrY, 'y')
      geneCumProb = self.geneCumProbPosY;
    elseif strcmp(posOrNeg, 'neg') && strcmp(xOrY, 'x')      
      geneCumProb = self.geneCumProbNegX;
    elseif strcmp(posOrNeg, 'neg') && strcmp(xOrY, 'y')
      geneCumProb = self.geneCumProbNegY;
    else
      error('Bug found');
    end
    % select axis
    if strcmp(xOrY, 'x')
      pc = self.coefXY(:,1);
    elseif strcmp(xOrY, 'y')
      pc = self.coefXY(:,2);
    else
      error('Bug found');
    end
    % threshold
    thr = (1-cutoff);    
    limit = self.solveLoadingLimit(geneCumProb, thr);        
    % select indices
    if strcmp(posOrNeg, 'pos')
      ind = find(pc >= limit);
    elseif strcmp(posOrNeg, 'neg')
      ind = find(pc <= -limit);
    else
      error('Bug found');
    end
    % merge with current indices
    ind = unique([ind(:)' self.geneSelIndices(:)']);
    self.genePm.setSelection(ind);  
  end
end

methods (Access = private)

end
  
end

