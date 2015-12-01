function varargout = pareto_plt_gui(varargin)
% PARETO_PLT_GUI MATLAB code for pareto_plt_gui.fig
%      PARETO_PLT_GUI, by itself, creates a new PARETO_PLT_GUI or raises the existing
%      singleton*.
%
%      H = PARETO_PLT_GUI returns the handle to a new PARETO_PLT_GUI or the handle to
%      the existing singleton*.
%
%      PARETO_PLT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARETO_PLT_GUI.M with the given input arguments.
%
%      PARETO_PLT_GUI('Property','Value',...) creates a new PARETO_PLT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pareto_plt_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pareto_plt_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pareto_plt_gui

% Last Modified by GUIDE v2.5 26-Apr-2015 23:16:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pareto_plt_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pareto_plt_gui_OutputFcn, ...
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

% --- Executes just before pareto_plt_gui is made visible.
function pareto_plt_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pareto_plt_gui (see VARARGIN)

% Choose default command line output for pareto_plt_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using pareto_plt_gui.
if strcmp(get(hObject,'Visible'),'off')
    if length(varargin)>0,main_data.d=varargin{1};else,main_data.d=[];end
    d=main_data.d;
    handles.pca_fig=figure;handles.pca_ax=gca;
    handles.gn_fig=figure;handels.gn_ax=gca;
    set(handles.pca_fig,'color','w');
    set(handles.gn_fig,'color','w');
    guidata(hObject,handles);
    %draw the pca plot window
    U=-1*leng
    [coeff,score,tsq,explained]=visualize_pca(nS,U,handles.pca_ax,handles.gn_ax,pcs)
    %variables to hold network centrality estimates
    %update gene data display on mouse move
    fh=@(varargin) mouse_mv_callback(handles.ax,pdata);
    set(handles.f,'windowbuttonmotionfcn',fh);
    main_data.pdata=pdata;main_data.oidx=oidx;
    main_data.glists_textbox=smp.glists_textbox;
    set(handles.pareto_plt_gui_root,'UserData',main_data);
end



% UIWAIT makes pareto_plt_gui wait for user response (see UIRESUME)
% uiwait(handles.pareto_plt_gui_root);


% --- Outputs from this function are returned to the command line.
function varargout = pareto_plt_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.pareto_plt_gui_root)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.pareto_plt_gui_root,'Name') '?'],...
                     ['Close ' get(handles.pareto_plt_gui_root,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.pareto_plt_gui_root)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});


% --- Executes on selection change in working_gene_list.
function working_gene_list_Callback(hObject, eventdata, handles)
% hObject    handle to working_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns working_gene_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from working_gene_list
main_data=get(handles.pareto_plt_gui_root,'UserData');
pdata=main_data.pdata;
glst_ids=get(handles.working_gene_list,'UserData');
idx=find(strcmpi(glst_ids{get(hObject,'Value')},pdata.lbls));
set(handles.gene_symbol,'String',pdata.lbls{idx});
set(handles.rank,'String',num2str(pdata.prank(idx)));
if pdata.fdr(idx)<=0.05
    set(handles.screen_fdr_textbox,'String',num2str(pdata.fdr(idx)),'ForegroundColor','r');
else
    set(handles.screen_fdr_textbox,'String',num2str(pdata.fdr(idx)),'ForegroundColor','k');
end
set(handles.nsh_textbox,'String',num2str(pdata.xdat(idx)));
set(handles.ntar_textbox,'String',num2str(pdata.ntar(idx)));
if isempty(pdata.tt)|pdata.tt(idx)==-1,
    set(handles.ttest_textbox,'String','NA');
    set(handles.fc_textbox,'String','NA');    
else
if pdata.tt(idx)<=0.05
    set(handles.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','r');
else
    set(handles.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','k');
end
set(handles.fc_textbox,'String',num2str(log2(pdata.fc(idx))));
end
set(handles.mlodz_textbox,'String',num2str(pdata.ydat(idx)));
set(handles.net_cent_textbox,'String',num2str(pdata.net_cent_comb(idx)));
set(handles.genetic_textbox,'String',num2str(pdata.net_cent_genetic(idx)));
set(handles.protein_textbox,'String',num2str(pdata.net_cent_protein(idx)));
set(handles.coloc_textbox,'String',num2str(pdata.net_cent_coloc(idx)));
set(handles.pathway_textbox,'String',num2str(pdata.net_cent_pathway(idx)));


% --- Executes during object creation, after setting all properties.
function working_gene_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to working_gene_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in select_genes.
function select_genes_Callback(hObject, eventdata, handles)
% hObject    handle to select_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~ishandle(handles.ax), return; end
jFigPeer = get(handle(gcf),'JavaFrame');
jWindow = jFigPeer.fHG1Client.getWindow;
title(handles.ax,sprintf(['Click or drag-click to add genes to the working list\n'...
                            'Ctrl+click or ctrl+drag-click to remove them\n'...
                            'Press enter key when done']),'FontSize',12)
main_data=get(handles.pareto_plt_gui_root,'UserData');
pdata=main_data.pdata;
figure(handles.f);
set(handles.f,'WindowStyle','modal');
%jf=get(get(handles.f,'JavaFrame'),'FigurePanelContainer');
%jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(1);
%jWindow.setEnabled(false);
h=gname_pareto(pdata.lbls(main_data.oidx),[],pdata);
set(handles.f,'WindowStyle','normal');
%jf.getComponent(0).getRootPane.getTopLevelAncestor.setAlwaysOnTop(0);
%jWindow.setEnabled(true);
if isempty(h),return;
else,title(handles.ax,'');end
%glst=get(handles.working_gene_list,'String');
%if strcmp('Add genes to the list',glst),glst={};end
for i=1:length(h)
    idx=min(find(strcmp(get(h(i),'String'),pdata.lbls)));
    gn=pdata.lbls{idx};
    glst_ids{i}=gn;
    if pdata.fdr(idx)<=0.05, gn=sprintf('<HTML><FONT color="%s">%s</font>', 'red', gn);end
    glst{i}=gn;
end
set(handles.working_gene_list,'String',glst,'UserData',glst_ids);




% --- Executes during object creation, after setting all properties.
function select_genes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in delete_selected_gene.
function delete_selected_gene_Callback(hObject, eventdata, handles)
% hObject    handle to delete_selected_gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glst=get(handles.working_gene_list,'String');
glst_ids=get(handles.working_gene_list,'UserData');
if isempty(glst_ids),return;end
n=length(glst_ids);
idx=get(handles.working_gene_list,'Value');
gn=glst_ids{idx};
glst_ids_new=glst_ids(1:idx-1);for i=idx+1:n,glst_ids_new{end+1}=glst_ids{i};end
if isempty(glst_ids_new)
    set(handles.working_gene_list,'String','Add genes to the list','UserData',[]);
else
    glst_new=glst(1:idx-1);for i=idx+1:n,glst_new{end+1}=glst{i};end
    set(handles.working_gene_list,'String',glst_new,'UserData',glst_ids_new);
end
h = findobj(handles.ax,'Type','text','Tag','gname');
i=1;while i<=length(h),if strcmp(get(h(i),'String'),gn),break;end;i=i+1;end
if ~isempty(h)&i>0&i<=length(h),delete(h(i));end
set(handles.working_gene_list,'Value',1);


% --- Executes on button press in add_gene.
function add_gene_Callback(hObject, eventdata, handles)
% hObject    handle to add_gene (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=get(handles.gene_to_find,'String');
main_data=get(handles.pareto_plt_gui_root,'UserData');
pdata=main_data.pdata;
idx=find(strcmpi(s,pdata.lbls));
if isempty(idx),set(handles.gene_to_find,'String','Gene not found!');return;end
glst=get(handles.working_gene_list,'String');
glst_ids=get(handles.working_gene_list,'UserData');
pidx=min(find(strcmpi(s,glst_ids)));
if ~isempty(pidx),
    figure(handles.f)
    moveptr(handle(handles.ax),'init');
    moveptr(handle(handles.ax),'move',pdata.xdat(idx),pdata.ydat(idx));
    set(handles.working_gene_list,'Value',pidx);
    set(handles.gene_symbol,'String',pdata.lbls{idx});
    set(handles.rank,'String',num2str(pdata.prank(idx)));
    if pdata.fdr(idx)<=0.05
        set(handles.fdr,'String',num2str(pdata.fdr(idx)),'ForegroundColor','r');
    else
        set(handles.fdr,'String',num2str(pdata.fdr(idx)),'ForegroundColor','k');
    end
    set(handles.nsh_textbox,'String',num2str(pdata.xdat(idx)));
    set(handles.ntar_textbox,'String',num2str(pdata.ntar(idx)));
    if isempty(pdata.tt)|pdata.tt(idx)==-1,
        set(handles.ttest_textbox,'String','NA');
        set(handles.fc,'String','NA');    
    else
        if pdata.tt(idx)<=0.05
            set(handles.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','r');
        else
            set(handles.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','k');
        end
        set(handles.fc_textbox,'String',num2str(log2(pdata.fc(idx))));
    end
    set(handles.mlodz_textbox,'String',num2str(pdata.ydat(idx)));
    set(handles.net_cent_textbox,'String',num2str(pdata.net_cent_comb(idx)));
    set(handles.genetic_textbox,'String',num2str(pdata.net_cent_genetic(idx)));
    set(handles.protein_textbox,'String',num2str(pdata.net_cent_protein(idx)));
    set(handles.coloc_textbox,'String',num2str(pdata.net_cent_coloc(idx)));
    set(handles.pathway_textbox,'String',num2str(pdata.net_cent_pathway(idx)));
    return;
elseif isempty(glst_ids),glst_ids={};glst={};end
figure(handles.f)
text(pdata.xdat(idx), pdata.ydat(idx),strjust(pdata.lbls{idx}), 'VerticalAlignment', 'baseline', 'Tag','gname','FontSize',12);
gn=pdata.lbls{idx};
glst_ids{end+1}=gn;
set(handles.gene_symbol,'String',gn);
set(handles.rank,'String',num2str(pdata.prank(idx)));
if pdata.fdr(idx)<=0.05
    set(handles.screen_fdr,'String',num2str(pdata.fdr(idx)),'ForegroundColor','r');
else
    set(handles.screen_fdr,'String',num2str(pdata.fdr(idx)),'ForegroundColor','k');
end
set(handles.nsh_textbox,'String',num2str(pdata.xdat(idx)));
set(handles.mlodz_textbox,'String',num2str(pdata.ydat(idx)));
set(handles.ntar_textbox,'String',num2str(pdata.ntar(idx)));
if isempty(pdata.tt)|pdata.tt(idx)==-1,
    set(handles.ttest_textbox,'String','NA');
    set(handles.fc_textbox,'String','NA');    
else
    if pdata.tt(idx)<=0.05
        set(handles.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','r');
    else
        set(handles.ttest_textbox,'String',num2str(pdata.tt(idx)),'ForegroundColor','k');
    end
    set(handles.fc_textbox,'String',num2str(log2(pdata.fc(idx))));
end
if pdata.fdr(idx)<=0.05, gn=sprintf('<HTML><FONT color="%s">%s</font>', 'red', gn);end
glst{end+1}=gn;
set(handles.working_gene_list,'String',glst,'UserData',glst_ids,'Value',length(glst));


function gene_to_find_Callback(hObject, eventdata, handles)
% hObject    handle to gene_to_find (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gene_to_find as text
%        str2double(get(hObject,'String')) returns contents of gene_to_find as a double


% --- Executes during object creation, after setting all properties.
function gene_to_find_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gene_to_find (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on gene_to_find and none of its controls.
function gene_to_find_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to gene_to_find (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in add_top_genes.
function add_top_genes_Callback(hObject, eventdata, handles)
% hObject    handle to add_top_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s1=get(handles.top_gene_cutoff,'String');
[cut1,ok]=str2num(s1);
if ~ok
    alert('title','Invalid FDR','String','Please enter a valid FDR, as a decimal. e.g. 0.05')
    return;
end
s2=get(handles.gene_rank_cutoff,'String');
[cut2,ok]=str2num(s2);
if ~ok
    alert('title','Invalid gene rank','String','Please enter a valid rank, as an integer. e.g. 100')
    return;
end
main_data=get(handles.pareto_plt_gui_root,'UserData');
pdata=main_data.pdata;
glst=get(handles.working_gene_list,'String');
glst_ids=get(handles.working_gene_list,'UserData');
if isempty(glst_ids),glst_ids={};glst={};end
idx=find(pdata.fdr<=cut1&pdata.prank<=cut2);
for i=1:length(idx)
    gn=pdata.lbls{idx(i)};
    if any(strcmpi(gn,glst_ids)),continue;end
    glst_ids{end+1}=gn;
    if pdata.fdr(idx(i))<=0.05, gn=sprintf('<HTML><FONT color="%s">%s</font>', 'red', gn);end
    glst{end+1}=gn;
    %figure(handles.f);
    %text(pdata.xdat(idx(i)), pdata.ydat(idx(i)),strjust(pdata.lbls{idx(i)}), 'VerticalAlignment', 'baseline', 'Tag','gname','FontSize',12);
end
set(handles.working_gene_list,'String',glst,'UserData',glst_ids,'Value',1);



function top_gene_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to top_gene_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of top_gene_cutoff as text
%        str2double(get(hObject,'String')) returns contents of top_gene_cutoff as a double


% --- Executes during object creation, after setting all properties.
function top_gene_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to top_gene_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_gene_list_to_file.
function save_gene_list_to_file_Callback(hObject, eventdata, handles)
% hObject    handle to save_gene_list_to_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.pareto_plt_gui_root,'UserData');
pdata=main_data.pdata;
glst=get(handles.working_gene_list,'String');
glst_ids=get(handles.working_gene_list,'UserData');
if isempty(glst_ids),return;end
[fname,pname, ~] = uiputfile({'*.txt','Text files (*.txt)';'*.tsv','Tab-separated value files (*.txv)'}, 'Save gene list as...');
if isequal(fname,0) || isequal(pname,0), return;end
f=fopen([pname fname],'w');
fprintf(f,['gene_symbol\trank\tscreen_p-value\tcollective_shRNA_activity_level\tactive_shRNA\tsequenced_shRNA\t'...
    'log2_expression_fold_change\texpression_t-test_p-value\n']);
for i=1:length(glst_ids)
    idx=min(find(strcmpi(glst_ids{i},pdata.lbls)));
    if isempty(idx),continue;end
    fprintf(f,'%s\t',pdata.lbls{idx});
    fprintf(f,'%d\t%0.5g\t%d\t%d\t%0.5g\t%0.5g\t%0.5g\n',...
        [pdata.prank(idx),pdata.fdr(idx),max(pdata.ydat(idx),0),pdata.xdat(idx),pdata.ntar(idx),log2(pdata.fc(idx)),pdata.tt(idx)]);
end
fclose(f);
alert('Title','Wrote gene list...','String',sprintf('Wrote current gene list to\n%s',fname));



% --- Executes on button press in gomap_pushbutton.
function gomap_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gomap_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.pareto_plt_gui_root,'UserData');
glst=get(handles.working_gene_list,'String');
glst_ids=get(handles.working_gene_list,'UserData');
main_glst=get(handles.main_glst,'String');%working gene list in main gui
main_glst_ids=get(handles.main_glst,'UserData');
glists=get(main_data.glists_textbox,'UserData');%gene lists textbox in main gui
glists_names=get(main_data.glists_textbox,'String');%names of gene lists in the main gui's gene lists textbox
lst_name=set_sample_id('title','Enter gene list ID:','string',sprintf(['Enter a name for the new gene list.']));
glists{end+1}=glst_ids;
if isempty(lst_name)||strcmp(lst_name,'Yes')
    glists_names{end+1}=['new_list' num2str(length(glists_names))];
else
    glists_names{end+1}=lst_name;
end
set(main_data.glists_textbox,'UserData',glists,'String',glists_names);
set(handles.main_glst,'String',glst_ids,'UserData',glst_ids,'Value',1);


% --- Executes during object creation, after setting all properties.
function net_cent_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to net_cent_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function gene_rank_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to gene_rank_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gene_rank_cutoff as text
%        str2double(get(hObject,'String')) returns contents of gene_rank_cutoff as a double


% --- Executes during object creation, after setting all properties.
function gene_rank_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gene_rank_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clear_gene_list_pushbutton.
function clear_gene_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clear_gene_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

glst=get(handles.working_gene_list,'String');
glst_ids=get(handles.working_gene_list,'UserData');
if isempty(glst_ids),return;end
h = findobj(handles.ax,'Type','text','Tag','gname');
for j=1:length(glst_ids)
    if isempty(h), break;end
    i=1;
    while i<=length(h)
        if ~ishandle(h(i))
            i=i+1;
            continue; 
        end
        if strcmp(get(h(i),'String'),glst_ids{j}),break;end;
        i=i+1;
    end
    if ~isempty(h)&i>0&i<=length(h),delete(h(i));end
end
set(handles.working_gene_list,'String',{'Add genes to the list'},'UserData',{});
set(handles.working_gene_list,'Value',1);



function pcax_Callback(hObject, eventdata, handles)
% hObject    handle to pcax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcax as text
%        str2double(get(hObject,'String')) returns contents of pcax as a double


% --- Executes during object creation, after setting all properties.
function pcax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pcay_Callback(hObject, eventdata, handles)
% hObject    handle to pcay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pcay as text
%        str2double(get(hObject,'String')) returns contents of pcay as a double


% --- Executes during object creation, after setting all properties.
function pcay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pcay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_pca_pushbutton.
function refresh_pca_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to refresh_pca_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cluster_cells_pushbutton.
function cluster_cells_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_cells_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in clustering_popup.
function clustering_popup_Callback(hObject, eventdata, handles)
% hObject    handle to clustering_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns clustering_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from clustering_popup


% --- Executes during object creation, after setting all properties.
function clustering_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clustering_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_trace_pushbutton.
function run_trace_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to run_trace_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in trace_popup.
function trace_popup_Callback(hObject, eventdata, handles)
% hObject    handle to trace_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trace_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trace_popup


% --- Executes during object creation, after setting all properties.
function trace_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trace_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitable4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function gene_annot_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gene_annot_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns ca
set(hObject,'Data',cell(4,1));


% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3


% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton26.
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton29.
function pushbutton29_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu14.
function popupmenu14_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu14


% --- Executes during object creation, after setting all properties.
function popupmenu14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15


% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitable3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
