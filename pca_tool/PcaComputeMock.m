classdef PcaComputeMock < PcaComputeBase
  %PCACOMPUTEMOCK Mock class to simulate real PcaCompute
  %   Detailed explanation goes here
  
  properties
    scoresXY
  end
  
%   properties (SetAccess = private, GetAccess = public)
%     n
%   end
  
  methods
    function self = PcaComputeMock()      
    end
    
    function computePca(self)
      n = size(self.d.nrmC, 2);
      self.coeff = rand(n,n);
      self.score = rand(n,n);
    end   
    
  end
  
end

