classdef PcaCompute < PcaComputeBase
  %PCACOMPUTE Computes PCA
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    function self = PcaCompute()
    end
    
    function computePca(self)
            if isempty(self.coeff) && ~isempty(self.d)
            [coeff,score,~,~,explained]=pca(self.d.nrmC');
            self.coeff = coeff';
            self.score = score';
            %self.explained=explained;
        end
    end
  end
  
end

