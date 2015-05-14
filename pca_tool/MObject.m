classdef MObject < handle
  %MOBJECT Signal and slot implemetation similar to Qt's
  %   See http://doc.qt.io/qt-5/signalsandslots.html
  %   In Matlab implementation, no need to explicitly define signals and
  %   slots.
  %   If an object want to emit a signal, it needs to be derived from 
  %   MObject. Just to receice signals (i.e. to contain slots) the
  %   derivation is not necessary.
  %
  
properties (GetAccess = public, SetAccess = private)
  signalToSlot
end

methods (Static)
  function connect(fromObj, signal, toObj, slot)
    fromObj.connectMe(signal, toObj, slot);
  end
end

methods
  function self = MObject()
    self.signalToSlot = containers.Map();
  end

  function connectMe(self, signal, toFuncH)
    self.signalToSlot(signal) = toFuncH;
  end
  
  function emit(self, signal, varargin)
    vals = self.signalToSlot(signal);
    N = length(varargin);
    for i=1:length(vals)
      toFuncH = vals(i);          
      copyArgs = cell(N,1);
      for j=1:N
        obj = varargin{j};
        if isa(obj, 'handle')
          copyArgs{j} = copy(obj);
        else
          copyArgs{j} = obj;
        end
      end
      % call the object method
      toFuncH(copyArgs{:});
    end    
  end  
  
%   function connectMe(self, signal, toObj, slot)
%     self.signalToSlot(signal) = {toObj, slot};
%   end
%   
%   function emit(self, signal, varargin)
%     vals = self.signalToSlot(signal);
%     N = length(varargin);
%     for i=1:size(vals,1)
%       toObj = vals{i,1};
%       slot = vals{i,2};          
%       copyArgs = cell(N,1);
%       for j=1:N
%         obj = varargin{j};
%         if isa(obj, 'handle')
%           copyArgs{j} = copy(obj);
%         else
%           copyArgs{j} = obj;
%         end
%       end
%       % call the object method
%       toObj.(slot)(copyArgs{:});
%     end    
%   end
end

end

