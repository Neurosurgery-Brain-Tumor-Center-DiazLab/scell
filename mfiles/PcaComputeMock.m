classdef PcaComputeMock < PcaComputeBase
  %PCACOMPUTEMOCK Mock class to simulate real PcaCompute
  %   Detailed explanation goes here
  
  properties (SetAccess = private, GetAccess = public)
    
  end
  
  methods
    function self = PcaComputeMock()      
    end
    
    function computePca(self)
      if ~isempty(self.d)
        geneCount = round(size(self.d.nrmC, 1)/10);
        sampleCount = size(self.d.nrmC, 2);
        span1 = 20*rand(1);
        span2 = 40*rand(1);
        self.coeff = span1*(rand(geneCount, geneCount) - 0.5);
        self.score = span2*(rand(sampleCount, sampleCount) -0.5);        
      end
    end   
    
    function computePcaUsingSamples(self, sampleIndices)  
      props = {'slbls', 'cidx', 'mapped', 'unmapped', 'ld_call', ...
        'turing', 'simpson', 'preseq', 'preseq_mar', 'lorenz', ...
        'lorenzh', 'sf', 'pareto', 'sidx'};
      d2 = self.d;
      for i = 1:length(props)
        p = props{i};
        tmp = d2.(p);
        d2.(p) = tmp(sampleIndices);
      end
      props = {'counts', 'cpm', 'nrmC'};
      for i = 1:length(props)
        p = props{i};
        tmp = d2.(p);
        d2.(p) = tmp(:, sampleIndices);
      end
      self.changeD(d2);
    end
    
  end
  
end

