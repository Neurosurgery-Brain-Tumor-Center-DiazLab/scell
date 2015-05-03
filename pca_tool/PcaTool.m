classdef PcaTool < GuiBase
  %PCATOOL PCA tool
  %   Detailed explanation goes here
  
  properties
    pcaFig
    scoreFig
    featureFig
  end
  
  properties (Constant)
    settingsFile = 'settings.mat'
    figHandles = {'pcaFig', 'scoreFig', 'featureFig'};
  end
  
  methods
    function self = PcaTool()
      self@GuiBase();
    end
    
    function show(self)
      if ishandle(self.pcaFig)
        set(self.pcaFig, 'Visible', 'on');
      else
        self.registerCallbacks();
        fig = pca_tool2();
        % save callbacks
        handles = guidata(fig);
        handles.n = {self};
        guidata(fig, handles);
        self.pcaFig = fig;
        self.createScatterPlots();
        self.loadSettings();
      end
    end
    
    %*** Callbacks from GUI
    function pareto_plt_gui_root_CloseRequestFcn(self)
      self.saveSettingsAndQuit();
    end
    
    %*** Other callbacks
    function saveSettingsAndQuit(self, varargin)
      self.saveSettings();
      % Delete all figs
      for nn = PcaTool.figHandles
        n = nn{1};
        if ishandle(self.(n))
          delete(self.(n));
        end
      end
    end
  end
  
  %*** Private implementation related stuff
  methods (Access = private)
    function createScatterPlots(self)
    % Creates scatter plot windows
      % pca scores
      fig = figure;
      set(fig, 'Name', 'Scores', 'CloseRequestFcn', ...
        @self.saveSettingsAndQuit);
      self.scoreFig = fig;
      % feature loadings
      fig = figure;
      set(fig, 'Name', 'Loadings', 'CloseRequestFcn', ...
        @self.saveSettingsAndQuit);      
      self.featureFig = fig;
    end
    
    function saveSettings(self)
    % Save settings, e.g, figure locations, etc., when closing the tool    
      settings = struct;      
      for i = 1:length(PcaTool.figHandles)
        n = PcaTool.figHandles{i};
        if ishandle(self.(n))
          settings.(n).outerPosition = get(self.(n), 'OuterPosition');
        end
      end
      save(PcaTool.settingsFile, 'settings')
    end
    
    function loadSettings(self)
    % Load settings when the tool is opened
    % For some reason, setting the GUI window position does not work
      if exist(PcaTool.settingsFile, 'file')
        tmp = load(PcaTool.settingsFile);
        settings = tmp.settings;
        fn = fieldnames(settings);
        for i = 1:length(fn)
          n = fn{i};
          if ishandle(self.(n))
            set(self.(n), 'OuterPosition', settings.(n).outerPosition);
          end
        end
      end
    end
  end
  
end

