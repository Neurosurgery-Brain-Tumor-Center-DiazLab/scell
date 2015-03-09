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

% Last Modified by GUIDE v2.5 23-Feb-2015 13:59:42

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
% uiwait(handles.gene_select_tool_root);
if strcmp(get(hObject,'Visible'),'off')
    main_data=get(handles.gene_select_tool_root,'UserData');
    if length(varargin)>0
        d=varargin{1};
        norm_tool_handle=varargin{2};
    else,return;end
    %test each gene for high variance (dispersion test)
    %test each gene for under-sampline (zero-inflation test)
    if ~isfield(d,'pnz')
        [pnz,zinf_fdr,iod,iod_fdr]=comp_gene_var_stats(d,true);
        d.pnz=pnz;d.zinf_fdr=zinf_fdr;
        d.iod=iod;d.iod_fdr=iod_fdr;
    end
    iod_fdr_cut=0.01;zinf_fdr_cut=0.01;pnz_cut=0.5;iod_cut=0.1;
    pl_iod=log(d.iod)/max(log(d.iod));
    idx1=find((d.iod_fdr<iod_fdr_cut&d.zinf_fdr>=zinf_fdr_cut)|(pl_iod>=(1-iod_cut)&d.pnz>=pnz_cut));
    idx2=setdiff(1:length(pl_iod),idx1);
    d.gidx=idx1;
    gidx=cell(length(pl_iod),1);
    set(handles.num_genes_text,'String',[num2str(length(idx1)) ' genes meet thresholds']);
    for i=1:length(idx1),gidx{idx1(i)}='Above thresholds';end
    for i=1:length(idx2),gidx{idx2(i)}='Below thresholds';end
    axes(handles.axes1);
    gscatter(d.pnz,pl_iod,gidx,'br','o',8);
    xlabel('Percent of cells expressing gene','FontSize',16);
    ylabel('Log index of dispersion percentile','FontSize',16)
    set(handles.axes1,'Xlim',[0 1],'Ylim',[0 1],'FontSize',16);
    l=findobj(gcf,'tag','legend');
    set(l,'location','SouthEast');
    main_data.d=d;
    main_data.norm_tool_handle=norm_tool_handle;
    set(handles.gene_select_tool_root,'UserData',main_data);
    [~,sidx]=sortrows(d.pnz,-1);
    [~,sidx2]=sortrows(d.iod(sidx),-1);
    sidx=sidx(sidx2);
    T=get(handles.gene_listbox,'Data');
    for i=1:length(d.gsymb)
        T{i,1}=d.gsymb{sidx(i)};
        T{i,2}=d.iod(sidx(i));
        T{i,3}=d.iod_fdr(sidx(i));
        T{i,4}=d.pnz(sidx(i));
        T{i,5}=d.zinf_fdr(sidx(i));
    end
    set(handles.gene_listbox,'Data',T);
end

% --- Outputs from this function are returned to the command line.
function varargout = gene_select_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
    set(handles.iod_edit,'String','0.1');
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
pl_iod=log(d.iod)/max(log(d.iod));
idx1=find((d.iod_fdr<iod_fdr_cut&d.zinf_fdr>=zinf_fdr_cut)|(pl_iod>=(1-iod_cut)&d.pnz>=pnz_cut));
idx2=setdiff(1:length(pl_iod),idx1);
set(handles.num_genes_text,'String',[num2str(length(idx1)) ' genes meet thresholds']);
d.gidx=idx1;
main_data.d=d;
set(handles.gene_select_tool_root,'UserData',main_data);
gidx=cell(length(pl_iod),1);
for i=1:length(idx1),gidx{idx1(i)}='Above thresholds';end
for i=1:length(idx2),gidx{idx2(i)}='Below thresholds';end
axes(handles.axes1);
gscatter(d.pnz,pl_iod,gidx,'br','o',8)
xlabel('Percent of cells expressing gene','FontSize',16);
ylabel('Log index of dispersion percentile','FontSize',16)
set(handles.axes1,'Xlim',[0 1],'Ylim',[0 1],'FontSize',16);
l=findobj(gcf,'tag','legend');
set(l,'location','SouthEast');


% --- Executes on button press in export_plot_button.
function export_plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to export_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.gene_select_tool_root,'UserData');
d=main_data.d;
%get thresholds from user
t=str2num(get(handles.iod_edit,'String'));
if isempty(t)
    iod_cut=0.1;
    set(handles.iod_edit,'String','0.1');
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
pl_iod=log(d.iod)/max(log(d.iod));
idx1=find(d.iod_fdr<iod_fdr_cut&d.zinf_fdr>=zinf_fdr_cut|pl_iod>=(1-iod_cut)|d.pnz>=pnz_cut);
idx2=setdiff(1:length(pl_iod),idx1);
gidx=cell(length(pl_iod),1);
for i=1:length(idx1),gidx{idx1(i)}='Above thresholds';end
for i=1:length(idx2),gidx{idx2(i)}='Below thresholds';end
figure
gscatter(d.pnz,pl_iod,gidx,'br','o',8)
xlabel('Percent of cells expressing gene','FontSize',16);
ylabel('Log index of dispersion percentile','FontSize',16)
set(gca,'Xlim',[0 1],'Ylim',[0 1],'FontSize',16);
l=findobj(gcf,'tag','legend');
set(l,'location','SouthEast');


% --- Executes on button press in export_gene_list_pushbutton.
function export_gene_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to export_gene_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T=get(handles.gene_listbox,'Data');
[fname pname]=uiputfile('*.txt','Save data as...','gene_list.txt');
try
    f=fopen(fullfile(pname,fname),'w');
catch me
    alert('String',['Error opening ' fullfile(pname,fname)]);
    return;
end
fprintf(f,'Gene\tIOD\tIOD_FDR\tPercent_expressing\tZero_inf_FDR\n');
for i=1:size(T,1)
    fprintf(f,'%s\t',T{i,1});
    fprintf(f,'%g\t',T{i,2});
    fprintf(f,'%g\t',T{i,3});
    fprintf(f,'%g\t',T{i,4});
    fprintf(f,'%g\n',T{i,5});
end
fclose(f);


% --- Executes on button press in select_genes.
function select_genes_Callback(hObject, eventdata, handles)
% hObject    handle to select_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.gene_select_tool_root,'UserData');
dnew=main_data.d;
norm_tool_data=get(main_data.norm_tool_handle,'UserData');
d=norm_tool_data.d;
t=str2num(get(handles.iod_edit,'String'));
if isempty(t)
    iod_cut=0.75;
    set(handles.iod_edit,'String','0.1');
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
pl_iod=log(dnew.iod)/max(log(dnew.iod));
d.gidx=find((dnew.iod_fdr<iod_fdr_cut&dnew.zinf_fdr>=zinf_fdr_cut)|(pl_iod>=(1-iod_cut)&dnew.pnz>=pnz_cut));
d.iod=dnew.iod;
d.pnz=dnew.pnz;
d.iod_fdr=dnew.iod_fdr;
d.zinf_fdr=dnew.zinf_fdr;
norm_tool_data.d=d;
set(main_data.norm_tool_handle,'UserData',norm_tool_data);
close(handles.gene_select_tool_root);
