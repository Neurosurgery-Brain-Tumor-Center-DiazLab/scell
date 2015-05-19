classdef UiState < handle
  %UISTATES Enumeration for different UI states
  %   One or more states can be active at the same time
  
  
properties (SetAccess = private, GetAccess = public)
  % toggle states
  genesSelected = false
  samplesSelected = false
  % buttons enabled based on UI state
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
  
  function updateGenesSelected(self, value)
    self.genesSelected = value;
    self.deleteGeneEnable = value;
    self.clearGeneListEnable = value;
    self.saveGeneListEnable = value;       
  end
  
  function updateSamplesSelected(self, value)
    self.samplesSelected = value;
    self.deleteSampleEnable = value;
    self.clearSampleListEnable = value;
    self.saveSampleListEnable = value; 
    self.refreshPcaUsingSamplesEnable = value;
  end
end


end

