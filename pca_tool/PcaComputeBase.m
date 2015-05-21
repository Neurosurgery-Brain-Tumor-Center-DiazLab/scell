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
  computePca(self)       
end

methods 
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

