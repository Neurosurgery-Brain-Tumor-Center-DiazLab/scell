classdef PcaComputeBase < handle & MObject
  %PCACOMPUTE Base class for computations needed by the PcaTool
  %   Detailed explanation goes here
  
properties (SetAccess = private, GetAccess = public)
  d
end

properties (SetAccess = protected, GetAccess = public)
  coeff
  score
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
end

end

