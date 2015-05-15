classdef PcaComputeMock < PcaComputeBase
  %PCACOMPUTEMOCK Mock class to simulate real PcaCompute
  %   Detailed explanation goes here
  
  properties (SetAccess = private, GetAccess = public)
    n
  end
  
  methods
    function self = PcaComputeMock()      
    end
    
    function computePca(self)
      if isempty(self.coeff) && ~isempty(self.d)
        n = size(self.d.nrmC, 2);
        self.coeff = rand(n,n);
        self.score = rand(n,n);
      end
    end   
    
  end
  
end

