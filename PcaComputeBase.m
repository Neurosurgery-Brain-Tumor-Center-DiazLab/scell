classdef PcaComputeBase < handle & MObject
  %PCACOMPUTE Base class for computations needed by the PcaTool
  %   Detailed explanation goes here
  
properties (SetAccess = protected, GetAccess = public)
  d
  coeff
  score
end

properties (Dependent)
  sampleCount
  geneCount
end

methods
  function self = PcaComputeBase()
  end
end

methods (Abstract)
  % Compute PCA given a matrix of normalized feature counts nrmC
  computePca(self, sampleIndices)       
  % PCA is recomputed using the given sample
  % indices. It's assumed that d exists
  computePcaUsingSamples(self, sampleIndices)
end

methods 
  function changeCoeff(self, cnew)
    self.coeff=cnew;
  end
  function changeScore(self, snew)
    self.score=snew;
  end
  function changeD(self, d)
    self.d = d;
    self.emit('d_changed');
  end
  
  function val = get.sampleCount(self)
    val = size(self.d.nrmC, 2);
  end
  
  function val = get.geneCount(self)
    val = size(self.d.nrmC, 1);
  end
end

end

