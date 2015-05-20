classdef ScatterSelectPM < MObject
  %SCATTERSELECTPM Presentation model for ScatterSelect
  %   Detailed explanation goes here
  
properties (GetAccess = public, SetAccess = private)
  data % Nx2 matrix of data points, each row = [x y]
  cluster % empty or Nx1 vector of cluster tags (integers >= 1), 
          % each corresponds to data row
  selData % selected data Mx2
  selIndices % incides to self.data corresponding to the selected data  
  currentIndex % index to data closest to the mouse button now, empty if
               % mouse is not inside the figure    
  clusterCount = 1 % number of clusters        
%   isInWindow = false % whether mouse pointer is inside the window
end

methods
  function self = ScatterSelectPM()
  end
    
  function updateData(self, data, cluster)
    self.data = data;
    self.clusterCount = 1;
    self.cluster = cluster;
    if ~isempty(cluster)
      tags = [];
      for i=1:length(cluster)
        if ~ismember(cluster(i), tags)
          tags(end+1) = cluster(i); %#ok<AGROW>
        end
      end
      self.clusterCount = length(tags);
    end
    self.emit('data_changed');
  end
  
  function updateCurrentIndex(self, value)
    self.currentIndex = value;
    self.emit('highlight_changed', value);    
  end  
  
  function selectCurrentIndex(self)
    if ~ismember(self.currentIndex, self.selIndices)
      self.selIndices(end+1) = self.currentIndex;
      self.selData = self.data(self.selIndices, :);
      self.emit('selection_changed');
    end    
  end
  
  function deselectIndex(self, index)
    self.selIndices(index) = [];
    self.selData(index,:) = [];
    self.emit('selection_changed');
  end
  
  function deselectLastIndex(self)
    if ~isempty(self.selIndices)
      self.deselectIndex(length(self.selIndices));
    end
  end
  
  function deselectAll(self)
    self.selIndices = [];
    self.selData = [];
    self.emit('selection_changed');
  end
  
end
  
end

