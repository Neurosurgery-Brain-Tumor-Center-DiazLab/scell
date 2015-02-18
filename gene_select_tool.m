function varargout = gene_select_tool(varargin)
% GENE_SELECT_TOOL MATLAB code for gene_select_tool.fig
%      GENE_SELECT_TOOL, by itself, creates a new GENE_SELECT_TOOL or raises the existing
%      singleton*.
%
%      H = GENE_SELECT_TOOL returns the handle to a new GENE_SELECT_TOOL or the handle to
%      the existing singleton*.
%
%      GENE_SELECT_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENE_SELECT_TOOL.M with the given input arguments.
%
%      GENE_SELECT_TOOL('Property','Value',...) creates a new GENE_SELECT_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gene_select_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gene_select_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gene_select_tool

% Last Modified by GUIDE v2.5 17-Feb-2015 17:58:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gene_select_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @gene_select_tool_OutputFcn, ...
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


% --- Executes just before gene_select_tool is made visible.
function gene_select_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gene_select_tool (see VARARGIN)

% Choose default command line output for gene_select_tool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gene_select_tool wait for user response (see UIRESUME)
% uiwait(handles.figure1);
if strcmp(get(hObject,'Visible'),'off')
    main_data=get(handles.gene_select_tool_root,'UserData');
    if length(varargin)>0,main_data.d=varargin{1};else,main_data.d=[];end
    set(handles.gene_select_tool_root,'UserData',main_data);
    [pnz,zinfp,zinf_fdr,iod,iodp,iod_fdr]=comp_gene_var_stats(d,true);
    d.pnz=pnz;d.zinfp=zinfp;d.zinf_fdr=zinf_fdr;
    d.iod=iod;d.iodp=iodp;d.iod_fdr=iod_fdr;
end


% --- Outputs from this function are returned to the command line.
function varargout = gene_select_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
main_data=get(handles.gene_select_tool_root,'UserData');
varargout{1} = main_data.d;



function iod_edit_Callback(hObject, eventdata, handles)
% hObject    handle to iod_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iod_edit as text
%        str2double(get(hObject,'String')) returns contents of iod_edit as a double


% --- Executes during object creation, after setting all properties.
function iod_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iod_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pnz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pnz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pnz_edit as text
%        str2double(get(hObject,'String')) returns contents of pnz_edit as a double


% --- Executes during object creation, after setting all properties.
function pnz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pnz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iod_fdr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to iod_fdr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iod_fdr_edit as text
%        str2double(get(hObject,'String')) returns contents of iod_fdr_edit as a double


% --- Executes during object creation, after setting all properties.
function iod_fdr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iod_fdr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zinf_fdr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zinf_fdr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zinf_fdr_edit as text
%        str2double(get(hObject,'String')) returns contents of zinf_fdr_edit as a double


% --- Executes during object creation, after setting all properties.
function zinf_fdr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zinf_fdr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_pushbutton.
function refresh_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.gene_select_tool_root,'UserData');
d=main_data.d;
%get thresholds from user
t=str2num(get(handles.iod_edit,'String'));
if isempty(t)
    iod_cut=0.75;
    set(handles.iod_edit,'String','0.75');
else
    iod_cut=t;
end
t=str2num(get(handles.pnz_edit,'String'));
if isempty(t)
    pnz_cut=0.5;
    set(handles.pnz_edit,'String','0.5');
else
    pnz_cut=t;
end
t=str2num(get(handles.iod_fdr_edit,'String'));
if isempty(t)
    iod_fdr_cut=0.01;
    set(handles.iod_fdr_edit,'String','0.01');
else
    iod_fdr_cut=t;
end
t=str2num(get(handles.zinf_fdr_edit,'String'));
if isempty(t)
    zinf_fdr_cut=0.01;
    set(handles.zinf_fdr_edit,'String','0.01');
else
    zinf_fdr_cut=t;
end
idx1=find(d.iod_fdr<0.01&d.zinf_fdr<0.01);
idx2=find(d.iod_fdr>0|d.zinf_fdr>0);
figure
scatter(d.pnz(idx1),log(d.iod(idx1))/log(max(d.iod)),'r')
hold
scatter(d.pnz(idx2),log(d.iod(idx2))/max(log(d.iod)),'b')
