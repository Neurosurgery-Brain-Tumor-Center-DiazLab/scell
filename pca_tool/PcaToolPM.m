classdef PcaToolPM < handle & MObject
  %PCATOOLPM Presentation model for the PcaTool
  %   Presentation model represent the logical contents of an UI.
  %   See http://martinfowler.com/eaaDev/PresentationModel.html
  
properties
  pcaxInd % pca x index
  pcayInd % pca y index
  clustering  % current ClusteringMethod
end
  
properties (SetAccess = private, GetAccess = public)
  compute % PcaComputeBase derived object
  coefXY % current pca axes coeff values 
  scoreXY % current pca axes score values
end

methods
  function self = PcaToolPM(computeObj)
    self@MObject();
    computeObj.connectMe('d_changed', @self.underlyingDataChanged);
    self.compute = computeObj;
    self.changeToSettings(self.defaultSettings);    
  end
  
  function val = defaultSettings(self)
    val.pcaxInd = 1;
    val.pcayInd = 2;
    val.clustering = ClusteringMethod.KMeans;
  end
  
  function val = getSettings(self)
    val.pcaxInd = self.pcaxInd;
    val.pcayInd = self.pcayInd;
    val.clustering = self.clustering;
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
      val = -1234;
    elseif strcmp(name, 'preseq_number')
      val = -1234;
    elseif strcmp(name, 'simpson_number')
      val = -1234;
    elseif strcmp(name, 'binomial_number')
      val = -1234;      
    else
      error('Unknown name %s', name);
    end
  end
  
  function changeToSettings(self, settings)
    self.pcaxInd = settings.pcaxInd;
    self.pcayInd = settings.pcayInd;
    self.clustering = settings.clustering;
  end
  
  function updateCurrentPca(self)
    self.compute.computePca();
    self.coefXY = [self.compute.coeff(:,self.pcaxInd) ...
                    self.compute.coeff(:,self.pcayInd)];
    self.scoreXY = [self.compute.score(:,self.pcaxInd) ...
                    self.compute.score(:,self.pcayInd)];                  
  end
  
  function underlyingDataChanged(self)
    self.emit('all_changed');
  end
  
end
  
end

