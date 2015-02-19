function varargout = norm_tool(varargin)
% NORM_TOOL MATLAB code for norm_tool.fig
%      NORM_TOOL, by itself, creates a new NORM_TOOL or raises the existing
%      singleton*.
%
%      H = NORM_TOOL returns the handle to a new NORM_TOOL or the handle to
%      the existing singleton*.
%
%      NORM_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NORM_TOOL.M with the given input arguments.
%
%      NORM_TOOL('Property','Value',...) creates a new NORM_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before norm_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to norm_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help norm_tool

% Last Modified by GUIDE v2.5 18-Feb-2015 13:02:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @norm_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @norm_tool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before norm_tool is made visible.
function norm_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to norm_tool (see VARARGIN)

% Choose default command line output for norm_tool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes norm_tool wait for user response (see UIRESUME)
% uiwait(handles.norm_tool_root);
if strcmp(get(hObject,'Visible'),'off')
    main_data=get(handles.norm_tool_root,'UserData');
    if length(varargin)>0,main_data.d=varargin{1};else,main_data.d=[];end
    set(handles.norm_tool_root,'UserData',main_data);
end


% --- Outputs from this function are returned to the command line.
function varargout = norm_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
main_data=get(handles.norm_tool_root,'UserData');
varargout{1} = main_data.d;


% --- Executes on button press in ercc_checkbox.
function ercc_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ercc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ercc_checkbox


% --- Executes on button press in bak_checkbox.
function bak_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to bak_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bak_checkbox


% --- Executes on button press in cyclin_checkbox.
function cyclin_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cyclin_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cyclin_checkbox


% --- Executes on button press in user_checkbox.
function user_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to user_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of user_checkbox


% --- Executes on button press in ercc_pushbutton.
function ercc_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ercc_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname pname]=uigetfile([],'Please load the ERCC gene identifiers.');
try
    D=importdata(fullfile(pname,fname));
catch me
    alert('String','File input error');
end
%find the ids in d and generate an index to them
main_data=get(handles.norm_tool_root,'UserData');
main_data.ercc_idx=[];
not_found={};
d=main_data.d;
for i=1:length(D)
    t=min(find(strcmp(D{i},d.gsymb)));
    if ~isempty(t), main_data.ercc_idx=[main_data.ercc_idx;t];
    else, not_found{end+1}=D{i};end
end
set(handles.norm_tool_root,'UserData',main_data);
if ~isempty(not_found)
    out_str=[num2str(length(not_found)),...
            sprintf(' genes could not be identified\n'),...
            sprintf('Genes not found:\n')];
    for u=1:length(not_found)-1,out_str=[out_str,sprintf('%s,',not_found{u})];end
    out_str=[out_str,sprintf('%s',not_found{end})];
    alert('String',out_str);
end 

% --- Executes on button press in bak_pushbutton.
function bak_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to bak_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cyclin_pushbutton.
function cyclin_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cyclin_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in user_pushbutton.
function user_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to user_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in select_genes_pushbutton.
function select_genes_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_genes_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.norm_tool_root,'UserData');
d=main_data.d;
dnew=gene_select_tool(d);
d.gidx=dnew.gidx;
d.iod=dnew.iod;
d.pnz=dnew.pnz;
d.iod_fdr=dnew.iod_fdr;
d.zinf_fdr=dnew.zinf_fdr;
main_data.d=d;
set(handles.norm_tool_root,'UserData',main_data);


% --- Executes on button press in norm_button.
function norm_button_Callback(hObject, eventdata, handles)
% hObject    handle to norm_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
