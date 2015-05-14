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
  signalToSlots
end

methods (Static)
  function connect(fromObj, signal, toObj, slot)
    fromObj.connectMe(signal, toObj, slot);
  end
end

methods
  function self = MObject()
    self.signalToSlots = containers.Map();
  end

  function connectMe(self, signal, toFuncH)
    if isKey(self.signalToSlots, signal)
      tmp = self.signalToSlots(signal);
      self.signalToSlots(signal) = [tmp {toFunch}];
    else      
      self.signalToSlots(signal) = {toFuncH};
    end
    
  end
  
  function emit(self, signal, varargin)
    if isKey(self.signalToSlots, signal)
      slots = self.signalToSlots(signal);
      N = length(varargin);
      for i=1:length(slots)
        toFuncH = slots{i};      
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
  end  
  
%   function connectMe(self, signal, toObj, slot)
%     self.signalToSlots(signal) = {toObj, slot};
%   end
%   
%   function emit(self, signal, varargin)
%     vals = self.signalToSlots(signal);
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

