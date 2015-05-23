classdef UiState < handle
  %UISTATES Enumeration for different UI states
  %   One or more states can be active at the same time
  
  
properties (SetAccess = private, GetAccess = public)
  % pca and clusterin
  pcxAxisEnable = false;
  pcyAxisEnable = false;
  refreshPcaEnable = false;
  clusteringMenuEnable = false;
  clusteringButtonEnable = false;  
  % gene and sample lists
  deleteGeneEnable = false
  clearGeneListEnable = false
  saveGeneListEnable = false
  deleteSampleEnable = false
  clearSampleListEnable = false
  saveSampleListEnable = false
  refreshPcaUsingSamplesEnable = false
end
  
methods
  function self = UiState()
  end
  
  function updateHasData(self, value)
    self.pcxAxisEnable = value;
    self.pcyAxisEnable = value;
    self.refreshPcaEnable = value;
    self.clusteringMenuEnable = value;
    self.clusteringButtonEnable = value;    
  end
  
  function updateGenesSelected(self, value)
    self.deleteGeneEnable = value;
    self.clearGeneListEnable = value;
    self.saveGeneListEnable = value;       
  end
  
  function updateSamplesSelected(self, value)
    self.deleteSampleEnable = value;
    self.clearSampleListEnable = value;
    self.saveSampleListEnable = value; 
    self.refreshPcaUsingSamplesEnable = value;
  end
end


end

