classdef PcaComputeMock2 < PcaComputeBase
  %PCACOMPUTEMOCK2 Mock class to simulate real PcaCompute with given PCA
  %   Detailed explanation goes here
  
  properties (SetAccess = private, GetAccess = public)
    
  end
  
  methods
    function self = PcaComputeMock2(coeff, score)      
      self.coeff = coeff;
      self.score = score;
    end
    
    function computePca(self)
      % does nothing
    end   
    
    function computePcaUsingSamples(self, sampleIndices)  
      props = {'slbls', 'cidx', 'mapped', 'unmapped', 'ld_call', ...
        'turing', 'simpson', 'preseq', 'preseq_mar', 'lorenz', ...
        'lorenzh', 'sf', 'pareto', 'sidx','ngns','ctype'};
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
      % coeff indices may go out of out with the current dataset
      [coeff,score]=pca(d2.nrmC');
      self.coeff = coeff;
      self.score = score;
      self.changeD(d2);
    end
    
  end
  
end

