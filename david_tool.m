function varargout = david_tool(varargin)
% DAVID_TOOL MATLAB code for david_tool.fig
%      DAVID_TOOL, by itself, creates a new DAVID_TOOL or raises the existing
%      singleton*.
%
%      H = DAVID_TOOL returns the handle to a new DAVID_TOOL or the handle to
%      the existing singleton*.
%
%      DAVID_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DAVID_TOOL.M with the given input arguments.
%
%      DAVID_TOOL('Property','Value',...) creates a new DAVID_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before david_tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to david_tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help david_tool

% Last Modified by GUIDE v2.5 10-Mar-2013 17:58:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @david_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @david_tool_OutputFcn, ...
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


% --- Executes just before david_tool is made visible.
function david_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to david_tool (see VARARGIN)

% % Choose default command line output for david_tool
% handles.output = hObject;
% % Update handles structure
% guidata(hObject, handles);
% 
% if length(varargin)>0,smp=varargin{1};else,smp=[];end
% main_data.id_type='ENTREZ gene ID';
% if strcmp(smp.genome{3},'hg19'),main_data.species='Homo sapien';
% elseif strcmp(smp.genome{3},'mm9'),main_data.species='Mus musculus';end
% main_data.fc=smp.fc;
% main_data.tt=smp.tt;
% main_data.pval=smp.pval;
% main_data.prank=smp.prank;
% main_data.nsh=smp.nsh;
% main_data.mlodz=smp.mlodz;
% main_data.gsymb=smp.gsymb;
% main_data.gid=smp.gid;
% main_data.eml=get(handles.edit2,'String');
% main_data.vis_dir='';
% main_data.pname=smp.pname;
% main_data.screen_pcut=smp.screen_pcut;
% main_data.gexp_pcut=smp.gexp_pcut;
% main_data.glists_textbox=smp.glists_textbox;
% main_data.main_glst=smp.main_glst;
% pareto_root_data=get(smp.pareto_gui_root_handle,'UserData');
% if ~pareto_root_data.java_loaded
%     h=waitbar(0,'loading java packaes. This only needs to be done once per session.');
%     d=dir('david_java_client/lib/*.jar');
%     wt=linspace(1/length(d),1,length(d));
%     for i=1:length(d)
%         s=fullfile(['david_java_client/lib/' d(i).name]);
%         waitbar(wt(i),h,{'Loading java package:',strrep(s,'_','\_')})
%         javaaddpath(s);   
%     end
%     import david.xsd.*;
%     import org.apache.axis2.AxisFault;
%     import sample.session.client.util.*;
%     import sample.session.service.xsd.*;
%     import sample.session.client.stub.*;
%     delete(h);
%     pareto_root_data.java_loaded=1;
%     set(smp.pareto_gui_root_handle,'UserData',pareto_root_data);
% end
% if isempty(pareto_root_data.GO),
%     wb=waitbar(0.5,'Loading gene ontology data. This needs to be done once per session.');
%     main_data.GO=geneont('file','gene_ontology.obo');
%     pareto_root_data.GO=main_data.GO;
%     set(smp.pareto_gui_root_handle,'UserData',pareto_root_data);
%     delete(wb);
% else
%     main_data.GO=pareto_root_data.GO;
% end
% %enter genes into the working list window
% for i=1:length(main_data.gsymb),gs{i}=main_data.gsymb{i};end
% set(handles.glist_textbox,'String',gs,'UserData',gs);
% set(handles.david_tool_root,'UserData',main_data);

% Choose default command line output for david_tool
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

smp=varargin{1};
main_data.id_type='Gene symbols';
main_data.species=smp.species; %set species and gsymb in caller
main_data.gsymb=smp.gsymb;
main_data.fc=[];main_data.tt=[];
main_data.pval=[];main_data.prank=[];
main_data.nsh=[];main_data.mlodz=[];
main_data.gid=[];

if strcmp(main_data.species,'hsa'),load('symb2entrez_hsa.mat');symb2ent=symb2ent_hsa;
else, load('symb2entrez_mmu.mat');symb2ent=symb2ent_mmu; end
kz=symb2ent.keys;
k=1;r=1;not_found={};
for i=1:length(main_data.gsymb) %look for each gene symbol's entrez id
    idx=find(strcmpi(main_data.gsymb{i},kz));
    if ~isempty(idx)
        en(k)=symb2ent(kz{min(idx)});
        k=k+1;
    else
        not_found{r}=main_data.gsymb{i};
        r=r+1;
    end                         
end
gene_data.gid=en;
if ~isempty(not_found)
    out_str=[num2str(length(not_found)),...
    sprintf(' genes could not be identified\n'),...
    sprintf('Genes not found:\n')];
    for u=1:length(not_found)-1,out_str=[out_str,sprintf('%s,',not_found{u})];end
    out_str=[out_str,sprintf('%s',not_found{end})];
    alert('String',out_str);
end 

main_data.eml=get(handles.edit2,'String');
main_data.vis_dir='';main_data.pname=[];
main_data.screen_pcut=[];main_data.gexp_pcut=[];
h=waitbar(0,'loading java packaes.');
d=dir('david_java_client/lib/*.jar');
wt=linspace(1/length(d),1,length(d));
for i=1:length(d)
    s=fullfile(['david_java_client/lib/' d(i).name]);
    waitbar(wt(i),h,{'Loading java package:',strrep(s,'_','\_')})
    javaaddpath(s);   
end
import david.xsd.*;
import org.apache.axis2.AxisFault;
import sample.session.client.util.*;
import sample.session.service.xsd.*;
import sample.session.client.stub.*;
delete(h);
wb=waitbar(0.5,'Loading gene ontology data.');
main_data.GO=geneont('file','gene_ontology.obo');
delete(wb);
set(handles.david_tool_root,'UserData',main_data);



% UIWAIT makes david_tool wait for user response (see UIRESUME)
% uiwait(handles.david_tool_root);


% --- Outputs from this function are returned to the command line.
function varargout = david_tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
main_data=get(handles.david_tool_root,'UserData');
contents = cellstr(get(hObject,'String'));
main_data.id_type=contents{get(hObject,'Value')};
set(handles.david_tool_root,'UserData',main_data);

% --- Executes during object creation, after setting all properties.
function glist_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to glist_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in query_david_pushbutton.
function query_david_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to query_david_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_data=get(handles.david_tool_root,'UserData');
if isempty(main_data.eml)
    alert('title','Enter an email in the edit box which is registered with DAVID','String','No valid email found.');
    return;
end
gene_data=get(handles.glist_textbox,'UserData');
c=query_david(gene_data.gid,main_data.eml);
x=pack_david_clusr_for_treemap(c,gene_data,main_data.GO);
gene_clusts={};
for i=1:length(x.children) %10 clusters at most, set the radio button string for each
    gene_clusts{i}=[];
    for j=1:length(x.children(i).children)
       D=textscan(x.children(i).children(j).data.gns,'%n','Delimiter',',');
       gene_clusts{i}=union(gene_clusts{i},D{1}); 
    end
    eval(['set(handles.c' num2str(i) '_button,''String'',x.children(' num2str(i) ').name,''Value'',1)']);
end
set(handles.all_button,'Value',1);
gene_data.gene_clusts=gene_clusts;
set(handles.glist_textbox,'UserData',gene_data);
main_data.gid=gene_data.gid;
main_data.gsymb=gene_data.gsymb;
if isfield(gene_data,'pval')
    main_data.pval=gene_data.pval;
    main_data.fc=gene_data.fc;
end
main_data.gene_disp_idx=1:length(gene_data.gsymb);
set(handles.david_tool_root,'UserData',main_data);
clear c;
pause(0.5);
if ~isdeployed
    outdir=uigetdir(main_data.pname,'Please select a directory where I can save your visualizations.');
else
    outdir=ctfroot;
end
if ~outdir, delete(h);return;end
main_data.vis_dir=outdir;
set(handles.david_tool_root,'UserData',main_data);
h=waitbar(0.25,'Writting javascript files');
js=make_treemap_json_from_david(x);
unzip(which('david_clustering.zip'),outdir);
f=fopen(fullfile(outdir,'david_clustering','data.js'),'w');
fprintf(f,'var json_data = %s;',js);
fclose(f);
waitbar(0.75,h,'Spawning web browser')
if ispc, dos(['start ' fullfile(outdir,'david_clustering','david_treemap.html') ' &']);
elseif ismac, unix(['open ' fullfile(outdir,'david_clustering','david_treemap.html') ' &']);
else unix(['firefox ' fullfile(outdir,'david_clustering','david_treemap.html') ' &']);end
%web(['file://' fullfile(outdir,'david_clustering','david_treemap.html')],'-browser')
if get(handles.to_file_radiobutton,'Value')
    waitbar(0.9,h,'Packaging web files for you to use later')
    [fname pname]=uiputfile(fullfile(outdir,'david_cluster_report.zip'));
    t=min(strfind(fname,'.'));
    if isempty(t),flnm=fname;else,flnm=fname(1:max(1,t-1));end
    if ispc,dos(['rename ' fullfile(outdir,'david_clustering') ' ' fullfile(outdir,flnm)])
    else, unix(['mv ' fullfile(outdir,'david_clustering') ' ' fullfile(pname,flnm)]),end
    if ~isequal(fname,0)&~isequal(pname,0),zip(fullfile(pname,fname),fullfile(pname,flnm));end
    if isdeployed,tfl=fullfile(ctfroot,'david_cluster_report.txt');
    else,tfl=fullfile(pwd,'david_cluster_report.txt');end
    f=fopen(tfl);
    g=fopen(fullfile(pname,[flnm '_goterm_list.tsv']),'w');
    D=textscan(f,'%s%n%n%s','Delimiter','\t');
    fprintf(g,'Term\tp-value\tPercent_annotated\tGenes\n');
    for i=1:length(D{1})
        fprintf(g,'%s\t',D{1}{i});
        fprintf(g,'%g\t',D{2}(i));
        fprintf(g,'%g\t',D{3}(i));
        fprintf(g,'%s\n',D{4}{i});
    end
    fclose(f);fclose(g);
end
delete(h)

% --- Executes on selection change in glist_textbox.
function glist_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to glist_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns glist_textbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from glist_textbox
gene_data=get(hObject,'UserData');%this holds entrez ids, pval, fold change, etc
gnz=get(hObject,'String');
idx=get(hObject,'Value');
if isfield(gene_data,'pval')&&~isempty(gene_data.pval)
    set(handles.ttest_textbox,'String',num2str(gene_data.pval(idx)));
    set(handles.exp_fc_textbox,'String',num2str(gene_data.fc(idx)));
else
    set(handles.ttest_textbox,'String','NA');
    set(handles.exp_fc_textbox,'String','NA');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
main_data=get(handles.david_tool_root,'UserData');
main_data.eml=get(hObject,'String');
set(handles.david_tool_root,'UserData',main_data);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in register_david_pushbutton.
function register_david_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to register_david_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('http://david.abcc.ncifcrf.gov/webservice/register.htm')

% --- Executes on button press in to_file_radiobutton.
function to_file_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to to_file_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of to_file_radiobutton


% --- Executes on button press in view_david_pushbutton.
function view_david_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to view_david_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.david_tool_root,'UserData');
if isempty(main_data.vis_dir),alert('Title','No visualizations found','String','You must run Query DAVID first');return;end
if ispc, dos(['start ' fullfile(main_data.vis_dir,'david_clustering','david_treemap.html')]);
elseif ismac, unix(['open ' fullfile(main_data.vis_dir,'david_clustering','david_treemap.html')]);
else unix(['firefox ' fullfile(main_data.vis_dir,'david_clustering','david_treemap.html')]);end


% --- Executes on button press in c1_button.
function c1_button_Callback(hObject, eventdata, handles)
% hObject    handle to c1_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c1_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{1})
    idxt=min(find(gene_data.gene_clusts{1}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c2_button.
function c2_button_Callback(hObject, eventdata, handles)
% hObject    handle to c2_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c2_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{2})
    idxt=min(find(gene_data.gene_clusts{2}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end
% --- Executes on button press in c3_button.
function c3_button_Callback(hObject, eventdata, handles)
% hObject    handle to c3_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c3_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{3})
    idxt=min(find(gene_data.gene_clusts{3}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c4_button.
function c4_button_Callback(hObject, eventdata, handles)
% hObject    handle to c4_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c4_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{4})
    idxt=min(find(gene_data.gene_clusts{4}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end


% --- Executes on button press in all_button.
function all_button_Callback(hObject, eventdata, handles)
% hObject    handle to all_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of all_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
if get(hObject,'Value')
    set(handles.glist_textbox,'String',main_data.gsymb);
    t.pval=main_data.pval;t.fc=main_data.fc;
    t.gsymb=main_data.gsymb;t.gid=main_data.gid;
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=1:length(main_data.gid);
    set(handles.david_tool_root,'UserData',main_data);
    for i=1:10,eval(['set(handles.c' num2str(i) '_button,''Value'',1)']);end
else
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'UserData',t,'Value',1);
    set(handles.glist_textbox,'String',{});
    main_data.gene_disp_idx=[];
    set(handles.david_tool_root,'UserData',main_data);
    for i=1:10,eval(['set(handles.c' num2str(i) '_button,''Value'',0)']);end
end

% --- Executes on button press in c5_button.
function c5_button_Callback(hObject, eventdata, handles)
% hObject    handle to c5_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c5_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{5})
    idxt=min(find(gene_data.gene_clusts{5}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c6_button.
function c6_button_Callback(hObject, eventdata, handles)
% hObject    handle to c6_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c6_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{6})
    idxt=min(find(gene_data.gene_clusts{6}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c7_button.
function c7_button_Callback(hObject, eventdata, handles)
% hObject    handle to c7_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c7_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{7})
    idxt=min(find(gene_data.gene_clusts{7}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c8_button.
function c8_button_Callback(hObject, eventdata, handles)
% hObject    handle to c8_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c8_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{8})
    idxt=min(find(gene_data.gene_clusts{8}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c9_button.
function c9_button_Callback(hObject, eventdata, handles)
% hObject    handle to c9_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c9_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{9})
    idxt=min(find(gene_data.gene_clusts{9}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in c10_button.
function c10_button_Callback(hObject, eventdata, handles)
% hObject    handle to c10_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c10_button
main_data=get(handles.david_tool_root,'UserData');
gene_data=get(handles.glist_textbox,'UserData');
if ~isfield(gene_data,'gene_clusts'),return;end
gl={};j=1;
for i=1:length(gene_data.gene_clusts{10})
    idxt=min(find(gene_data.gene_clusts{10}(i)==main_data.gid));
    if ~isempty(idxt),idx(j)=idxt;j=j+1;end
end
get(hObject,'Value')
if get(hObject,'Value')
    cur_idx=main_data.gene_disp_idx;
    cur_idx=union(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
else
    cur_idx=main_data.gene_disp_idx;
    cur_idx=setdiff(cur_idx,idx);
    t=[];t.gsymb=main_data.gsymb(cur_idx);t.gid=main_data.gid(cur_idx);
    if ~isempty(main_data.pval)
        t.pval=main_data.pval(cur_idx);t.fc=main_data.fc(cur_idx);
    end
    t.gene_clusts=gene_data.gene_clusts;
    set(handles.glist_textbox,'String',main_data.gsymb(cur_idx));
    set(handles.glist_textbox,'UserData',t,'Value',1);
    main_data.gene_disp_idx=cur_idx;
    set(handles.david_tool_root,'UserData',main_data);
    set(handles.all_button,'Value',0);
end

% --- Executes on button press in delete_gene_pushbutton.
function delete_gene_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to delete_gene_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx=get(handles.glist_textbox,'Value');
glist_old=get(handles.glist_textbox,'String');
glist_data_old=get(handles.glist_textbox,'UserData');
if isempty(idx)
    set(handles.glist_textbox,'String',{},'UserData',[],'Value',1);
else
    j=1;
    for i=1:idx-1
        glist_data.gid(j)=glist_data_old.gid(i);
        glist_data.gsymb(j)=glist_data_old.gsymb(i);
        if isfield(glist_data_old,'pval')&&~isempty(glist_data_old.pval)
            glist_data.pval(j)=glist_data_old.pval(i);
            glist_data.fc(j)=glist_data_old.fc(i);
            glist_data.prank(j)=glist_data_old.prank(i);
        end
        glist{j}=glist_old{i};
        j=j+1;
    end
    for i=idx+1:length(glist_old)
        glist_data.gid(j)=glist_data_old.gid(i);
        glist_data.gsymb(j)=glist_data_old.gsymb(i);
        if isfield(glist_data_old,'pval')&&~isempty(glist_data_old.pval)
            glist_data.pval(j)=glist_data_old.pval(i);
            glist_data.fc(j)=glist_data_old.fc(i);
            glist_data.prank(j)=glist_data_old.prank(i);
        end
        glist{j}=glist_old{i};
        j=j+1;
    end
    set(handles.glist_textbox,'String',glist,'UserData',glist_data,'Value',max(1,idx-1));
end




function find_gene_editbox_Callback(hObject, eventdata, handles)
% hObject    handle to find_gene_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of find_gene_editbox as text
%        str2double(get(hObject,'String')) returns contents of find_gene_editbox as a double


% --- Executes during object creation, after setting all properties.
function find_gene_editbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to find_gene_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_glist_pushbutton.
function save_glist_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_glist_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gns=get(handles.glist_textbox,'String');
gene_data=get(handles.glist_textbox,'UserData');
main_data=get(handles.david_tool_root,'UserData');
if ~isfield(main_data,'last_dir')
    [fname pname]=uiputfile('*.*','Save gene list to file...');
else
    [fname pname]=uiputfile([main_data.last_dir,'*.*'],'Save gene list to file...');
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.david_tool_root,'UserData',main_data);
    f=fopen(fullfile(pname,fname),'w');
end
for i=1:length(gns)
    fprintf(f,'%s\t',gns{i});
    if ~isempty(gene_data.pval)
        fprintf(f,'%g\t',gene_data.fc(i));
        fprintf(f,'%g\t',gene_data.pval(i));
    fprintf(f,'\n');
    end
end
alert('String',['Wrote ' fname '...']);

% --- Executes on button press in find_gene_button.
function find_gene_button_Callback(hObject, eventdata, handles)
% hObject    handle to find_gene_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=get(handles.find_gene_editbox,'String');
if isempty(s), return;end
gene_data=get(handles.glist_textbox,'UserData');
idx=min(find(strcmpi(s,gene_data.gsymb)));
if isempty(idx),set(handles.find_gene_editbox,'String','Gene not found!');return;end
if ~isempty(idx)
    set(handles.glist_textbox,'Value',idx);
    if isempty(gene_data.pval)
        set(handles.ttest_textbox,'String','NA');
        set(handles.exp_fc_textbox,'String','NA');
    else
        set(handles.ttest_textbox,'String',num2str(gene_data.pval(idx)));
        set(handles.exp_fc_textbox,'String',num2str(gene_data.fc(idx)));
    end
end

% --- Executes during object creation, after setting all properties.
function query_david_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to query_david_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in delete_list_pushbutton.
function delete_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to delete_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
idx=get(handles.gene_lists,'Value');
if length(glists_names)==1
    set(handles.gene_lists,'UserData',{},'String',{});
    set(handles.glist_textbox,'UserData',[],'String',{});
    return;
else
    j=1;
    for i=1:idx-1
        gl_data{j}=glists_data{i};
        gl_names{j}=glists_names{i};
        j=j+1;
    end
    for i=idx+1:length(glists_data)
        gl_names{j}=glists_names{i};
        gl_data{j}=glists_data{i}
        j=j+1;
    end
    set(handles.gene_lists,'UserData',gl_data,'String',gl_names,'Value',max(1,idx-1));
    t=gl_data{max(1,idx-1)};
    set(handles.glist_textbox,'UserData',t,'String',t.gsymb);
end
% --- Executes on button press in export_gene_list_pushbutton.
function export_gene_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to export_gene_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sample_data=get(handles.sample_list,'UserData');%get screen data
main_data=get(handles.david_tool_root,'UserData');%get global data
glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
if isempty(glists_data),return;end
idx=get(handles.gene_lists,'Value');
gnz=glists_data{idx};
if ~isfield(main_data,'last_dir')
    [fname pname]=uiputfile(glists_names{idx},'Select a file to write to...');
else
    main_data.last_dir
    [fname pname]=uiputfile('*.*','Select a file to write to...',fullfile(main_data.last_dir,[strrep(strrep(glists_names{idx},' ','_'),',',''),'.txt']));
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.david_tool_root,'UserData',main_data);
end
f=fopen(fullfile(pname,fname),'w');
fprintf(f,'gene\tscreen_rank\tscreen_pvalue\tcollective_hairpin_activity\tnumber_of_active_hairpins');
if ~isempty(sample_data.fc)
    fprintf(f,'\ttreatment_over_control_fold_change\tttest_pvalue');
end
if ~isempty(sample_data.net_cent_comb)
    fprintf(f,'\tnetwork_centrality');
end
fprintf(f,'\n');
for i=1:length(gnz)
    pidx=min(find(strcmpi(gnz{i},sample_data.gsymb)));
    fprintf(f,'%s\t',gnz{i});
    fprintf(f,'%u\t',sample_data.prank(pidx));
    fprintf(f,'%g\t',sample_data.pval(pidx));
    fprintf(f,'%g\t',sample_data.mlodz(pidx));
    fprintf(f,'%u',sample_data.nsh(pidx));
    if ~isempty(sample_data.fc)
        fprintf(f,'\t%g\t',sample_data.fc(pidx));
        fprintf(f,'%g',sample_data.tt(pidx));
    end
    if ~isempty(sample_data.net_cent_comb),fprintf(f,'\t%g',sample_data.net_cent_comb(pidx));end
    fprintf(f,'\n');
end

% --- Executes on button press in new_list_pushbutton.
function new_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to new_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
gene_data=get(handles.glist_textbox,'UserData');
glists_data{end+1}=gene_data;
lst_name=set_sample_id('title','Enter gene list ID:','string',sprintf(['Enter a name for the gene list.']));
if isempty(lst_name),glists_names{end+1}=['new_list' num2str(length(glists_names))];
else,glists_names{end+1}=lst_name;end
set(handles.gene_lists,'UserData',glists_data,'String',glists_names);


% --- Executes on button press in import_gene_list_pushbutton.
function import_gene_list_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to import_gene_list_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

main_data=get(handles.david_tool_root,'UserData');%get global data
glists_data=get(handles.gene_lists,'UserData');%get gene lists data
glists_names=get(handles.gene_lists,'String');%get the name of all the gene lists
idx=get(handles.gene_lists,'Value');
list_info=choose_file_type('title','Select data type','String','Select gene ID type and organism');
main_data.id_type=list_info.sample_type;
main_data.species=list_info.genome;
isent=strcmp(main_data.id_type,'ENTREZ IDs');
ge=list_info.ge;
if isent, s='%n';else, s='%s';end
if ge,s=[s '%n%n'];end
s=[s '%*[^\n]'];
%get the filename
if ~isfield(main_data,'last_dir')
    [fname pname]=uigetfile('*.*','Select the gene list...');
else
    [fname pname]=uigetfile([main_data.last_dir,'*.*'],'Select the gene list...');
end
if ~isstr(fname),return;
else
    main_data.last_dir=pname;
    set(handles.david_tool_root,'UserData',main_data);
    f=fopen(fullfile(pname,fname));
end
D=textscan(f,s,'Delimiter','\t');
D{1}=deblank(D{1});
if ge
    fc=D{2};%gene_data.fc=D{2};
    pval=D{3};%gene_data.pval=D{3};
    [~,sidx]=sort(abs(fc),'descend');
    prank=zeros(size(sidx));
    prank(sidx)=1:length(sidx);%gene_data.prank=sidx;
end
if isent
    gene_data.gid=D{1};gs={};
    if strcmp(main_data.species,'Other')
        for i=1:length(D{1}),gs{i}=num2str(D{1}(i));end
    else
        if strcmp(main_data.species,'Homo sapiens')
            load('entrez2gsymb_hsa.mat');ent2symb=ent2symb_hsa;
        end
        if strcmp(main_data.species,'Mus musculus')
            load('entrez2gsymb_mmu.mat');ent2symb=ent2symb_mmu;
        end
        for i=1:length(D{1})
            if ent2symb.isKey(D{1}(i)),gs{i}=ent2symb(D{1}(i));
            else, gs{i}=num2str(D{1}(i)),end
        end
    end
    gene_data.gsymb=gs;
    if exist('fc','var')
        gene_data.fc=fc;
        gene_data.pval=pval;
        gene_data.prank=prank;
    end
else
    isref=strcmp(main_data.id_type,'RefSeq accession');
    ishum=strcmp(main_data.species,'Homo sapiens');
    if isref
        if ishum,load('refseq2entrez_hsa.mat');symb2ent=ref2ent_hsa;
        else,load('refseq2entrez_mmu.mat');symb2ent=ref2ent_mmu;
        end
    else
        if ishum,load('symb2entrez_hsa.mat');symb2ent=symb2ent_hsa;
        else,load('symb2entrez_mmu.mat');symb2ent=symb2ent_mmu;
        end
    end
    kz=symb2ent.keys;
    k=1;r=1;not_found={};
    for i=1:length(D{1}) %look for each gene symbol's entrez id
        idx=find(strcmpi(D{1}{i},kz));
        if ~isempty(idx)
            gs{k}=D{1}{i};
            en(k)=symb2ent(kz{min(idx)});
            if exist('fc','var')
                gene_data.fc(k)=fc(i);
                gene_data.pval(k)=pval(i);
                gene_data.prank(k)=prank(i);
            end
            k=k+1;
        else
            not_found{r}=D{1}{i};r=r+1;
        end                         
    end
    gene_data.gid=en;
    gene_data.gsymb=gs;
        if ~isempty(not_found)
            out_str=[num2str(length(not_found)),...
                sprintf(' genes could not be identified\n'),...
                sprintf('try using unique identifiers like ENTREZ ids. Genes not found:\n')];
            for u=1:length(not_found)-1,out_str=[out_str,sprintf('%s,',not_found{u})];end
            out_str=[out_str,sprintf('%s',not_found{end})];
            alert('String',out_str);
        end 
end
set(handles.glist_textbox,'String',gs,'UserData',gene_data,'Value',1);
lst_name=[];
lst_name=set_sample_id('title','Enter gene list ID:','string',sprintf(['Enter a name for the gene list\n(' fname ')']));
if isempty(lst_name),glists_names{end+1}=['new_list' num2str(length(glists_names))];
else,glists_names{end+1}=lst_name;end
glists_data{end+1}=gene_data;
set(handles.gene_lists,'UserData',glists_data,'String',glists_names,'Value',length(glists_names));
set(handles.david_tool_root,'UserData',main_data);


% --- Executes on selection change in gene_lists.
function gene_lists_Callback(hObject, eventdata, handles)
% hObject    handle to gene_lists (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns gene_lists contents as cell array
%        contents{get(hObject,'Value')} returns selected item from gene_lists
glists_data=get(hObject,'UserData');%get gene lists data
glists_names=get(hObject,'String');%get the name of all the gene lists
idx=get(hObject,'Value');
if ~isempty(idx)
    gene_data=glists_data{idx};
    set(handles.glist_textbox,'String',gene_data.gsymb,'UserData',glists_data{idx},'Value',1);
    main_data=get(handles.david_tool_root,'UserData');
    if isfield(gene_data,'pval')
        main_data.pval=gene_data.pval;
        main_data.fc=gene_data.fc;
        main_data.prank=gene_data.prank;
    else
        main_data.pval=[];
        main_data.fc=[];
        main_data.prank=[];
    end
    main_data.gsymb=gene_data.gsymb;
    main_data.gid=gene_data.gid;
    set(handles.david_tool_root,'UserData',main_data);
    for i=1:10
        eval(['set(handles.c' num2str(i) '_button,''String'',''Cluster ' num2str(i) ''',''Value'',0)']);
    end
    set(handles.all_button,'Value',0);
else
    set(handles.glist_textbox,'String',{},'UserData',{},'Value',1);
    main_data.pval=[];main_data.fc=[];main_data.prank=[];
    main_data.gsymb=[];main_data.gid=[];
    set(handles.david_tool_root,'UserData',main_data);
    for i=1:10
        eval(['set(handles.c' num2str(i) '_button,''String'',''Cluster ' num2str(i) ''',''Value'',0)']);
    end
    set(handles.all_button,'Value',0);
end


% --- Executes during object creation, after setting all properties.
function gene_lists_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gene_lists (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
