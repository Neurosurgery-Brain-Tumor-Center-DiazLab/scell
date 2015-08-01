classdef ScatterSelectPM < MObject
  %SCATTERSELECTPM Presentation model for ScatterSelect
  %   Detailed explanation goes here
  
properties (GetAccess = public, SetAccess = private)
  title = 'Unnamed' % figure title
  data % Nx2 matrix of data points, each row = [x y]
  cluster % empty or Nx1 vector of cluster tags (integers >= 1), 
          % each corresponds to data row
  selData % selected data Mx2
  selIndices % incides to self.data corresponding to the selected data  
  currentIndex % index to data closest to the mouse button now, empty if
               % mouse is not inside the figure    
  clusterCount = 1 % number of clusters        
%   isInWindow = false % whether mouse pointer is inside the window
  selectionMode = SelectionMode.Pointwise % selection mode
  lassoColor = 'magenta'
  lassoPoints % Nx2 matrix of (x,y) lasso polygon points, the last point
              % is assumed to be connected to the first point
end

methods
  function self = ScatterSelectPM()
  end
    
  
  function set.title(self, value)
    self.title = value;
    self.emit('changed');
  end

  function updateData(self, data, cluster)
    self.data = data;
    self.cluster = cluster;
    if isempty(cluster)
      self.clusterCount = 1;
    else
      self.clusterCount = length(unique(cluster));
    end
    self.selData = self.data(self.selIndices,:);
    self.emit('data_changed');
  end
  
  function addLassoPoint(self, xy)
    self.lassoPoints(end+1,:) = xy(:)';
    self.emit('changed');
  end
  
  function undoLassoPoint(self, isAll)
    if ~isempty(self.lassoPoints)
      if isAll
        self.lassoPoints = [];
      else
        self.lassoPoints(end,:) = [];        
      end
      self.emit('changed');
    end
  end
  
  function tf = selectLassoedPoints(self)
  % Select points inside the lassoed polygon using the winding number
  % algorithm http://geomalgorithms.com/a03-_inclusion.html
    tf = false;
    if ~isempty(self.lassoPoints) && size(self.lassoPoints,1) > 2 ...
        && ~isempty(self.data)
      tf = true;
      N = size(self.data,1);
      V = [self.lassoPoints; self.lassoPoints(1,:)];
      insideIndices = [];
      for i = 1:N
        P = self.data(i,:);
        if wn_PnPoly(P, V) ~= 0
          insideIndices(end+1) = i; %#ok<AGROW>
        end
      end
      newIndices = unique([self.selIndices insideIndices]);
      self.setSelection(newIndices);
    end
  end
  
  function updateSelectionMode(self, mode)
    if self.selectionMode ~= mode
      self.selectionMode = mode;
      self.emit('changed');
    end
  end
  
  function updateCurrentIndex(self, value)
    self.currentIndex = value;
    self.emit('highlight_changed', value);    
  end  
  
  function selectCurrentIndex(self)
    self.selectIndex(self.currentIndex);
  end
  
  function selectIndex(self, index)
    if ~ismember(index, self.selIndices)
      self.selIndices(end+1) = index;
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
  
  function setSelection(self, indices)
    self.selIndices = indices;
    self.selData = self.data(indices,:);
    self.emit('selection_changed');
  end
  
  function deselectAll(self)
    self.selIndices = [];
    self.selData = [];
    self.emit('selection_changed');
  end
  
end
  
end

%*** 3rd party code adapted for matlab
% http://geomalgorithms.com/a03-_inclusion.html

% Copyright 2000 softSurfer, 2012 Dan Sunday
% This code may be freely used, distributed and modified for any purpose
% providing that this copyright notice is included with it.
% SoftSurfer makes no warranty for this code, and cannot be held
% liable for any real or imagined damage resulting from its use.
% Users of this code must verify correctness for their application.

function val = isLeft(P0, P1, P2)
% // isLeft(): tests if a point is Left|On|Right of an infinite line.
% //    Input:  three points P0, P1, and P2
% //    Return: >0 for P2 left of the line through P0 and P1
% //            =0 for P2  on the line
% //            <0 for P2  right of the line
% //    See: Algorithm 1 "Area of Triangles and Polygons"
% inline int
% isLeft( Point P0, Point P1, Point P2 )
% {
%     return ( (P1.x - P0.x) * (P2.y - P0.y)
%             - (P2.x -  P0.x) * (P1.y - P0.y) );
% }
  
  % Matlab: P0 = [x y], etc.

  val = ( (P1(1) - P0(1)) * (P2(2) - P0(2)) ...
              - (P2(1) - P0(1)) * (P1(2) - P0(2)) );
end

function wn = wn_PnPoly(P, V)
% // wn_PnPoly(): winding number test for a point in a polygon
% //      Input:   P = a point,
% //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
% //      Return:  wn = the winding number (=0 only when P is outside)
% int
% wn_PnPoly( Point P, Point* V, int n )
% {
%     int    wn = 0;    // the  winding number counter
% 
%     // loop through all edges of the polygon
%     for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
%         if (V[i].y <= P.y) {          // start y <= P.y
%             if (V[i+1].y  > P.y)      // an upward crossing
%                  if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
%                      ++wn;            // have  a valid up intersect
%         }
%         else {                        // start y > P.y (no test needed)
%             if (V[i+1].y  <= P.y)     // a downward crossing
%                  if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
%                      --wn;            // have  a valid down intersect
%         }
%     }
%     return wn;
% }

  % Matlab: P = [x y], V = Nx2 matrix of vertices [x y] 
  wn = 0;
  n = size(V,1)-1;
  % loop through all edges of the polygon
  for i = 1:n
    if V(i,2) <= P(2)
      if V(i+1,2) > P(2)
        if isLeft(V(i,:), V(i+1,:), P) > 0
          wn = wn + 1;
        end
      end
    else
      if V(i+1,2) <= P(2)
        if isLeft(V(i,:), V(i+1,:), P) < 0
          wn = wn -1;
        end
      end
    end
  end

end
