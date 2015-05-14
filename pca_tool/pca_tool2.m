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

% Last Modified by GUIDE v2.5 14-May-2015 15:50:22

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


% --- Executes on button press in findGeneButtonH.
function findGeneButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to findGeneButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in refreshPcaButtonH.
function refreshPcaButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to refreshPcaButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function pcaxEditH_Callback(hObject, eventdata, handles)
% hObject    handle to pcaxEditH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcaxEditH as text
%        str2double(get(hObject,'String')) returns contents of pcaxEditH as a double



function pcayEditH_Callback(hObject, eventdata, handles)
% hObject    handle to pcayEditH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcayEditH as text
%        str2double(get(hObject,'String')) returns contents of pcayEditH as a double


% --- Executes on selection change in clusteringPopupH.
function clusteringPopupH_Callback(hObject, eventdata, handles)
% hObject    handle to clusteringPopupH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clusteringPopupH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clusteringPopupH


% --- Executes on button press in clusterCellsButtonH.
function clusterCellsButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to clusterCellsButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tracePopupH.
function tracePopupH_Callback(hObject, eventdata, handles)
% hObject    handle to tracePopupH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tracePopupH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tracePopupH


% --- Executes on button press in runTraceButtonH.
function runTraceButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to runTraceButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in geneListboxH.
function geneListboxH_Callback(hObject, eventdata, handles)
% hObject    handle to geneListboxH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns geneListboxH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from geneListboxH



function geneSymbolEditH_Callback(hObject, eventdata, handles)
% hObject    handle to geneSymbolEditH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of geneSymbolEditH as text
%        str2double(get(hObject,'String')) returns contents of geneSymbolEditH as a double


% --- Executes on button press in saveGeneListButtonH.
function saveGeneListButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to saveGeneListButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in addTopGenesButtonH.
function addTopGenesButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to addTopGenesButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in selectGenesButtonH.
function selectGenesButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to selectGenesButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in deleteGeneButtonH.
function deleteGeneButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to deleteGeneButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearGeneListButtonH.
function clearGeneListButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to clearGeneListButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ontologyButtonH.
function ontologyButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to ontologyButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function cutoffEditH_Callback(hObject, eventdata, handles)
% hObject    handle to cutoffEditH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoffEditH as text
%        str2double(get(hObject,'String')) returns contents of cutoffEditH as a double


% --- Executes on selection change in pc1PopupH.
function pc1PopupH_Callback(hObject, eventdata, handles)
% hObject    handle to pc1PopupH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pc1PopupH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pc1PopupH


% --- Executes on selection change in posPopupH.
function posPopupH_Callback(hObject, eventdata, handles)
% hObject    handle to posPopupH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns posPopupH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from posPopupH


% --- Executes on selection change in sampleListboxH.
function sampleListboxH_Callback(hObject, eventdata, handles)
% hObject    handle to sampleListboxH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns sampleListboxH contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sampleListboxH


% --- Executes on button press in findSampleButtonH.
function findSampleButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to findSampleButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function sampleSymbolEditH_Callback(hObject, eventdata, handles)
% hObject    handle to sampleSymbolEditH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampleSymbolEditH as text
%        str2double(get(hObject,'String')) returns contents of sampleSymbolEditH as a double


% --- Executes on button press in saveSampleListButtonH.
function saveSampleListButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to saveSampleListButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in selectSamplesButtonH.
function selectSamplesButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to selectSamplesButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in deletSampleButtonH.
function deletSampleButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to deletSampleButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearSampleListButtonH.
function clearSampleListButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to clearSampleListButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in refreshPcaUsingSamplesButtonH.
function refreshPcaUsingSamplesButtonH_Callback(hObject, eventdata, handles)
% hObject    handle to refreshPcaUsingSamplesButtonH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
