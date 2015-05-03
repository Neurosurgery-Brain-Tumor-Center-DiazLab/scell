function varargout = pca_tool2(varargin)
%PCA_TOOL2 M-file for pca_tool2.fig
%      PCA_TOOL2, by itself, creates a new PCA_TOOL2 or raises the existing
%      singleton*.
%
%      H = PCA_TOOL2 returns the handle to a new PCA_TOOL2 or the handle to
%      the existing singleton*.
%
%      PCA_TOOL2('Property','Value',...) creates a new PCA_TOOL2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to pca_tool2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PCA_TOOL2('CALLBACK') and PCA_TOOL2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PCA_TOOL2.M with the given input
%      arguments.
% 
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pca_tool2

% Last Modified by GUIDE v2.5 03-May-2015 15:30:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pca_tool2_OpeningFcn, ...
                   'gui_OutputFcn',  @pca_tool2_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
    % perform callbacks on the controller object too    
    if ~isempty(varargin) && ~isempty(varargin{4}) 
      handles = varargin{4};
      funcName = varargin{1};
      % there can be several callback objects
      for i=1:length(handles.n)
        handles.n{i}.guiCallback(funcName, varargin{2:end});        
      end
    end
end
% End initialization code - DO NOT EDIT


% --- Executes just before pca_tool2 is made visible.
function pca_tool2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for pca_tool2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pca_tool2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pca_tool2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when user attempts to close pareto_plt_gui_root.
function pareto_plt_gui_root_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to pareto_plt_gui_root (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
