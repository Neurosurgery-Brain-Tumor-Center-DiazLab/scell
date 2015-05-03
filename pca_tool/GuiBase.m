classdef GuiBase < dynamicprops
  %GUIBASE Base class fot Matlab GUIs
  %   Provides flexible handles and callbacks
  
  properties (SetAccess = private, GetAccess = public)
    callbackMap
  end
  
  methods
    function self = GuiBase()
    end
    
    function registerCallbacks(self)
    % Register all callback functions, i.e., those public methods that have
    % specified strings, e.g., '_Fcn' or '_Callback', in their name
      self.callbackMap = containers.Map();            
      m = metaclass(self);
      for i=1:length(m.MethodList)
        method = m.MethodList(i);
        if (~isempty(strfind(method.Name, '_Callback')) ...
           || ~isempty(strfind(method.Name, 'Fcn')))...
            && strcmp(method.Access, 'public')
          % number of input parameters
          parIn = length(method.InputNames)-1;
          % save
          self.callbackMap(method.Name) = {...
            @(varargin)self.(method.Name)(varargin{:}),...
            parIn};
        end
      end      
    end
    
    function guiCallback(self, name, varargin)
    % Reroute callbacks from GUI
    % varargin{1} == hObject, i.e., event source handle
    % varargin{2} == eventdata
    % varargin{3} == handles
      if isKey(self.callbackMap, name)
        val = self.callbackMap(name);
        fh = val{1};
        parIn = val{2};
        fh(varargin{1:parIn});        
      else
        disp(['Callback ''' name ''' not found']);
      end
    end 
    
    function mapUiControls(self, parent)
    % Map matching user interface controls under 'parent' to properties 
    % of this object
    % 'parent' is usually a figure handle    
      children = findall(parent);
      for i=1:length(children)
        tag = get(children(i), 'Tag');
        if ~isempty(tag)
          if ~isprop(self, tag);
            addprop(self, tag);            
          end      
          self.(tag) = children(i);
        end
      end      
    end
  end
end

