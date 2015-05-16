classdef PcaComputeMock < PcaComputeBase
  %PCACOMPUTEMOCK Mock class to simulate real PcaCompute
  %   Detailed explanation goes here
  
  properties (SetAccess = private, GetAccess = public)
    
  end
  
  methods
    function self = PcaComputeMock()      
    end
    
    function computePca(self)
      if isempty(self.coeff) && ~isempty(self.d)
        geneCount = round(size(self.d.nrmC, 1)/10);
        sampleCount = size(self.d.nrmC, 2);
        self.coeff = rand(geneCount, geneCount);
        self.score = rand(sampleCount, sampleCount);
      end
    end   
    
  end
  
end

