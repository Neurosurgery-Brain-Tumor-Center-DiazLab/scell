classdef PcaToolPM < handle & MObject
  %PCATOOLPM Presentation model for the PcaTool
  %   Presentation model represent the logical contents of an UI.
  %   See http://martinfowler.com/eaaDev/PresentationModel.html
  
properties (SetAccess = private, GetAccess = public)
  pcaxInd % pca x index
  pcayInd % pca y index
  compute % PcaComputeBase derived object
  coefXY % current pca axes coeff values 
  scoreXY % current pca axes score values
end

methods
  function self = PcaToolPM(computeObj)
    self@MObject();
    computeObj.connectMe('d_changed', @self.underlyingDataChanged);
    self.compute = computeObj;
    self.changeToSettings(self.defaultSettings);    
  end
  
  function val = defaultSettings(self)
    val.pcaxInd = 1;
    val.pcayInd = 2;    
  end
  
  function val = getSettings(self)
    val.pcaxInd = self.pcaxInd;
    val.pcayInd = self.pcayInd;
  end
  
  function changeToSettings(self, settings)
    self.pcaxInd = settings.pcaxInd;
    self.pcayInd = settings.pcayInd;
  end
  
  function updateCurrentPca(self)
%     self.compute.computePca(
  end
  
  function underlyingDataChanged(self)
    self.emit('all_changed');
  end
  
end
  
end

