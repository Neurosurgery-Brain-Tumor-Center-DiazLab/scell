classdef PcaCompute < PcaComputeBase
  %PCACOMPUTE Computes PCA
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function self = PcaCompute()
    end
    
    function computePca(self)
      if ~isempty(self.d)
      % call the pca.m from here and set properties as in PcaComputeMock
      end
    end
  end
  
end

